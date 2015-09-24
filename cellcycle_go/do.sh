## -- created on Aug 11, 2015 --
## -- cell cycle related GOs are : --
## -- GO:0007049 cell cycle  --> GO:0000075 cell cycle checkpoint
## --                        --> GO:0051726 regulation of cell cycle
## --                        --> GO:0031570 DNA integrity checkpoint

perl ~/Dropbox/perl_scripts/gene_ontology/go_children.pl \
	-i GO:0007049 \
	-g /Users/wchen/databases/gene_ontology/june2015/go.obo \
	> cell_cycle.go

perl ~/Dropbox/perl_scripts/gene_ontology/go_children.pl \
	-i GO:0000075 \
	-g /Users/wchen/databases/gene_ontology/june2015/go.obo \
	> cell_cycle_checkpoint.go

perl ~/Dropbox/perl_scripts/gene_ontology/go_children.pl \
	-i GO:0051726 \
	-g /Users/wchen/databases/gene_ontology/june2015/go.obo \
	> regulation_of_cell_cycle.go

perl ~/Dropbox/perl_scripts/gene_ontology/go_children.pl \
	-i GO:0031570 \
	-g /Users/wchen/databases/gene_ontology/june2015/go.obo \
	> DNA_integrity_checkpoint.go
	
perl ~/Dropbox/perl_scripts/gene_ontology/go_children.pl \
	-i GO:0006281 \
	-g /Users/wchen/databases/gene_ontology/june2015/go.obo \
	> DNA_repair.go
	 