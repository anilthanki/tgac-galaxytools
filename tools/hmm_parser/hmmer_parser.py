"""
Simple parser to convert a HMMER search -tblout into hcluster_sg style output
"""
import argparse
import math
import sqlite3
import tempfile
import re


def create_tables(conn):
    cur = conn.cursor()
    cur.execute('''CREATE TABLE hmmer (
        id INTEGER PRIMARY KEY,
        sequence_id VARCHAR NOT NULL,
        hmm_id VARCHAR NOT NULL,
        evalue REAL NOT NULL)''')
    conn.commit()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='in-file', type=argparse.FileType('rt'), required=True, help='Path to input file')
    parser.add_argument('-o', metavar='out-file', type=argparse.FileType('wt'), required=True, help='Path to output file')
    
    options = parser.parse_args()

    db_file = tempfile.NamedTemporaryFile(suffix=".sqlite")

    conn = sqlite3.connect(db_file.name)
    conn.execute('PRAGMA foreign_keys = ON')

    create_tables(conn)

    cur = conn.cursor()

    for line in options.i:
        line = line.rstrip()
        if line.startswith("#"):
            continue
        else:
            line_cols = re.split(r'\s+', line)

        sequence_id = line_cols[0]
        hmm_id = line_cols[2]
        evalue = float(line_cols[4])

        cur.execute('SELECT evalue from hmmer WHERE sequence_id = ?', (sequence_id, ))

        row = cur.fetchone()

        if not row:
            cur.execute('INSERT INTO hmmer (sequence_id, hmm_id, evalue) VALUES (?, ?, ?)',
                    (sequence_id, hmm_id, evalue))
        else:
            if row[0] < evalue:
                cur.execute('UPDATE hmmer set evalue = ? WHERE sequence_id = ?',
                    (evalue, sequence_id))

        conn.commit()

    options.i.close()

    query = 'SELECT sequence_id, hmm_id, evalue FROM hmmer ORDER BY id'
    cur.execute(query)

    cluster_dict = {}

    while True:
        rows = cur.fetchone()

        if not rows:
            break
        else:
            if rows[1] in cluster_dict:
                cluster_dict[rows[1]].append(rows[0])
            else:
                cluster_dict[rows[1]] = [rows[0]]

    s = ","
    
    for x in cluster_dict:
        options.o.write("%s\t%d\t%s\n" % (x, len(cluster_dict[x]), s.join(cluster_dict[x])))

    conn.close()
    db_file.close()
    options.o.close()


if __name__ == "__main__":
    main()
