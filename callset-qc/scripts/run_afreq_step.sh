
OUT_DIR=/storage/nicholema/callset_qc_outs/afreq_step
SAMPLES=/gymreklab-tscc/mousavi/analysis/1000genomes/samples/all_samples.txt

GANGSTR_OUTS=/gymreklab-tscc/mousavi/results/1000genomes/outs/


mkdir -p $OUT_DIR

COMMA_SEP_STATS=""
SUPER_POPS=$(cat $SAMPLES | cut -f3 | sort | uniq)
for SPOP in $SUPER_POPS
do
	SPOP_DIR=$OUT_DIR/$SPOP
	mkdir -p $SPOP_DIR
	POPS=$(cat $SAMPLES | grep $SPOP | cut -f2 | sort | uniq)
	for POP in $POPS
	do
		POP_DIR=$SPOP_DIR/$POP
		mkdir -p $POP_DIR
		STAT_FILE=$POP_DIR/$POP.statstr.tab
		if [ -f "$STAT_FILE" ]; then
		    echo "$SPOP $POP already copied"
		else 
		    echo "Copying $SPOP $POP"
		    cp $GANGSTR_OUTS/$SPOP/$POP/merged/${POP}_merged.stats.tab $STAT_FILE
		fi
		COMMA_SEP_STATS="$COMMA_SEP_STATS,$STAT_FILE"		
	done
done
COMMA_SEP_STATS="${COMMA_SEP_STATS:1}"


cd ..
./generate_afreqs.sh > $OUT_DIR/known_afreqs.txt
./compile_afreqs.sh $COMMA_SEP_STATS $OUT_DIR/known_afreqs.txt $OUT_DIR > $OUT_DIR/called_afreqs.txt
