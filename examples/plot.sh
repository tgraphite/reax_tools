input=$1
output=${input/dot/svg}
if [[ -f $output ]]; then
    rm $output
fi
dot -Tsvg $input > $output