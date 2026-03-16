echo "# Parameters"
echo
echo "\`\`\`"
ganon -h
echo "\`\`\`"
modes=( "build" "build-custom" "update" "classify" "reassign" "report" "table" )
for m in "${modes[@]}"
do
    echo
    echo "<details>"
    echo "  <summary>ganon ${m}</summary>"
    echo
    echo "\`\`\`"
    ganon ${m} -h
    echo "\`\`\`"
    echo
    echo "</details>"
done
