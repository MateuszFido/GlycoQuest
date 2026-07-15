## apache
sudo apt-get install apache2 
sudo apt-get install libapache-session-perl 
## other
sudo apt-get install bioperl 
sudo apt-get install libberkeleydb-perl
sudo apt-get install libcgi-fast-perl
sudo apt-get install libcgi-session-perl
sudo apt-get install libdata-dumper-concise-perl
sudo apt-get install libdata-dumper-simple-perl
sudo apt-get install libhtml-template-perl
sudo apt-get install libfile-copy-recursive-perl
sudo apt-get install libgd-graph-perl
sudo apt-get install libio-compress-bzip2-perl
sudo apt-get install libio-compress-perl
sudo apt-get install libtemplate-perl
sudo apt-get install libxml-treebuilder-perl
sudo apt-get install libxml-writer-perl
sudo apt-get install libmldbm-perl
sudo apt-get install libstatistics-descriptive-perl
sudo apt-get install libcgi-formbuilder-perl
sudo apt-get install libmail-sender-perl
sudo apt-get install build-essential
sudo apt-get install tofrodos
sudo apt-get install subversion

## create softlink for dos2unix
cd /usr/bin
sudo ln -s fromdos dos2unix
cd ~

## via cpan
sudo cpan CGI::FastTemplate
sudo cpan HTML:PageIndex
sudo cpan Math:Random
#sudo cpan Mail::Sender
#sudo cpan CGI::Ajax
#sudo cpan Email::Valid