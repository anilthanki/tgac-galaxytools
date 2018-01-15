"""
Simple parser to convert a BLAST 12-column or 24-column tabular output into a
3-column tabular input for hcluster_hg (id1, id2, weight):
"""
import argparse
import math
import sqlite3
import tempfile
from time import gmtime, strftime
print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))


BATCH_SIZE = 2000


def create_tables(conn):
    cur = conn.cursor()
    cur.execute('''CREATE TABLE alignment (
        id INTEGER PRIMARY KEY,
        sequence1_id VARCHAR NOT NULL,
        sequence2_id VARCHAR NOT NULL,
        weight DOUBLE NOT NULL,
        score DOUBLE NOT NULL,
        sequence1_self_score DOUBLE NOT NULL,
        sequence2_self_score DOUBLE NOT NULL,
        blast_score_ratio DOUBLE NOT NULL,
        reciprocal BOOLEAN NOT NULL
        )''')
    conn.commit()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='in-file', type=argparse.FileType('rt'), required=True, help='Path to input file')
    parser.add_argument('-o', metavar='out-file', type=argparse.FileType('wt'), required=True, help='Path to output file')
    parser.add_argument('-r', action='store_true', default=False,
                        dest='reciprocal',
                        help='Returns only reciprocal results')
    parser.add_argument('-e', default=0,
                        dest='ensembl',
                        help='Returns results which are either reciprocal or having BLAST Score Ratio more than 0.33')
    options = parser.parse_args()

    db_file = tempfile.NamedTemporaryFile(suffix=".sqlite")

    conn = sqlite3.connect(db_file.name)
    conn.execute('PRAGMA foreign_keys = ON')

    create_tables(conn)

    cur = conn.cursor()

    i = 0
    for line in options.i:
        line = line.rstrip()
        line_cols = line.split('\t')
        sequence1_id = line_cols[0]
        sequence2_id = line_cols[1]

        # Ignore self-matching hits
        # if sequence1_id == sequence2_id:
            # continue

        i = i + 1
        evalue = float(line_cols[10])


        # Convert evalue to an integer weight with max 100
        if evalue != 0.0:
            weight = min(100, round(math.log10(evalue) / -2.0))
        else:
            # If the evalue is 0, leave weight at 100
            weight = 100

        score = line_cols[13]
        cur.execute('INSERT INTO alignment (id, sequence1_id, sequence2_id, weight, score, blast_score_ratio, sequence1_self_score, sequence2_self_score, reciprocal) VALUES (?, ?, ?, ?, ?, ?, ? ,? ,?)',
                    (i, sequence1_id, sequence2_id, weight, score, 0, 0, 0, 0))

        # Commit transaction at every BATCH_SIZE rows to save memory
        if i % BATCH_SIZE == 0:
            conn.commit()

    # Commit final transaction
    conn.commit()
    options.i.close()

    # Delete alternative alignments keeping only one with max weight
    cur.execute('DELETE FROM alignment WHERE id NOT IN (SELECT id FROM alignment GROUP BY sequence1_id, sequence2_id HAVING weight=max(weight))')
    conn.commit()

    # Find sequence1_self_score
    cur.execute('UPDATE alignment SET sequence1_self_score = (SELECT a.score FROM alignment a WHERE a.sequence1_id = alignment.sequence1_id and a.sequence1_id = a.sequence2_id)')
    conn.commit()

    # Find sequence2_self_score
    cur.execute('UPDATE alignment SET sequence2_self_score = (SELECT a.score FROM alignment a WHERE a.sequence1_id = alignment.sequence2_id and a.sequence1_id = a.sequence2_id)')
    conn.commit()

     # Delete self matching alignments
    cur.execute('DELETE FROM alignment WHERE sequence1_id = sequence2_id')
    conn.commit()

    # Calculate BLAST Score ratio
    cur.execute('UPDATE alignment SET blast_score_ratio = score/MAX(sequence1_self_score, sequence2_self_score)')
    conn.commit()

    # Updating Reciprocal results
    cur.execute('UPDATE alignment SET reciprocal = 1 where id in (SELECT a1.id FROM alignment a1, alignment a2 WHERE a1.sequence1_id = a2.sequence2_id AND a1.sequence2_id = a2.sequence1_id ORDER BY a1.id)')
    conn.commit()

    if options.ensembl:
        if options.ensembl is "1":
            query = 'SELECT sequence1_id, sequence2_id, weight FROM alignment WHERE reciprocal = 1 OR blast_score_ratio > 0.33'
        elif options.ensembl is "2":
            query = 'SELECT sequence1_id, sequence2_id, weight FROM alignment WHERE blast_score_ratio > 0.11'
    elif options.reciprocal:
        query = 'SELECT sequence1_id, sequence2_id, weight FROM alignment WHERE reciprocal = 1'
    else:
        query = 'SELECT sequence1_id, sequence2_id, weight FROM alignment ORDER BY id'

    cur.execute(query)
    while True:
        rows = cur.fetchmany(BATCH_SIZE)
        if not rows:
            break
        for row in rows:
            options.o.write("%s\t%s\t%d\n" % row)

    conn.close()
    db_file.close()
    options.o.close()
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))


if __name__ == "__main__":
    main()
