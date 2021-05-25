detex project.tex | tr -d '[:punct:]' | tr -d '[:digit:]' | tr ' ' '\n' | sort | uniq -c | sort -rn > words
