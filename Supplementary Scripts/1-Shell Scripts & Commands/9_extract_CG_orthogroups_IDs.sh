awk -F'\t' '
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i == "Homo_sapiens") col = i
  }
}
NR>1 && $col ~ /\S/ {
  print $1
}
' Orthogroups.tsv > HCG_orthogroups_IDs.txt

