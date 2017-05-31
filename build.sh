make clean
rm -f very-simple-theory-replication.zip
rm -f /tmp/very-simple-theory-replication.zip

working_directory=$PWD

cd ..
Rscript -e "Rcpp::compileAttributes('very-simple-theory-replication/actyR')"
Rscript -e "devtools::document('very-simple-theory-replication/actyR')"
Rscript -e "devtools::test('very-simple-theory-replication/actyR')"

cd very-simple-theory-replication
./render.rb README.md > README.html
cd ..

zip -r /tmp/very-simple-theory-replication.zip very-simple-theory-replication \
	--exclude "*.so" "*.o" "*build.sh" "*/\.*" "*.RData" "*.Rproj" "*.tex" "*.yml" "*.pdf" "*.log" "*README.md" "*render.rb" "*.css"
cd $working_directory

rm -rf /tmp/very-simple-theory-replication
mv -f /tmp/very-simple-theory-replication.zip .
