INSTALLDIR=/home/xqxp/xquest/V2_1_1/xquest

chmod -R 755 $INSTALLDIR/bin

## Configuration for the Webserver
echo "Webserver configuration ..."
chmod 755 $INSTALLDIR/cgi/changeheader.pl
chmod -v 755 $INSTALLDIR/cgi/*.cgi
echo "updating cgi headers"
cd $INSTALLDIR/cgi
./changeheader.pl
cd $INSTALLDIR
find . -type f -exec dos2unix {} \;

## chmod the logfolder
chmod -R 777 $INSTALLDIR/logs

## cp sylesheet to the webfolder
read -p "Do you want to copy the css stylesheet to a location that can be accessed by the webserver?(yes/no)" exportvar
if [ "$exportvar" = "yes" ]; then
read -p "Please specify location (/var/www/)" csslocation
if [ $csslocation ] 
then
echo "Will copy style.css to $csslocation"
else
csslocation=/var/www/
fi
echo "cp $INSTALLDIR/styles/style.css $csslocation/style.css"
sudo cp $INSTALLDIR/styles/style.css $csslocation/style.css
fi

## add /xquest/bin directory to PATH
echo "Please add the path to the xquest/bin directory manually to your PATH if you install xquest the first time!"