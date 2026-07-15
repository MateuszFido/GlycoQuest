./bin/mergexml.pl -list resultdirectories_fullpath -resdir $1 -v
cd $1
cp merged_xquest.xml xquest.xml
../bin/annotatexml.pl -in xquest.xml -out annotated_xquest.xml -native -v
rm xquest.xml
../bin/xprophet.pl
../bin/xprophet.pl -in annotated_xquest.xml -out xquest.xml
cd ..
cp -R $1 /IMSB/ra/mewing/html/xquestresults/
chmod -R 777 /IMSB/ra/mewing/html/xquestresults/$1
