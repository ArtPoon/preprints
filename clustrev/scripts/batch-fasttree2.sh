for f in `ls ../data/*.*_TRUE.fas`; do fasttree2 -nt -gtr < $f > $f.tree; done

