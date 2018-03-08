CREATE TABLE transcript (
        transcript_id VARCHAR PRIMARY KEY,
        protein_id VARCHAR UNIQUE,
        protein_sequence VARCHAR,
        gene_id VARCHAR NOT NULL REFERENCES gene(gene_id));
CREATE TABLE meta (
        version VARCHAR);
CREATE TABLE gene_family_member (
        gene_family_id INTEGER NOT NULL REFERENCES gene_family(gene_family_id),
        protein_id VARCHAR KEY NOT NULL REFERENCES transcript(protein_id),
        protein_alignment VARCHAR NOT NULL,
        PRIMARY KEY (gene_family_id, protein_id));
CREATE TABLE gene_family (
        gene_family_id INTEGER PRIMARY KEY,
        gene_tree VARCHAR NOT NULL);
CREATE TABLE gene (
        gene_id VARCHAR PRIMARY KEY,
        symbol VARCHAR,
        species VARCHAR NOT NULL,
        gene_json VARCHAR NOT NULL);
CREATE VIEW transcript_species as
        SELECT transcript_id, species 
        FROM transcript JOIN gene
        ON gene.gene_id = transcript.gene_id;
