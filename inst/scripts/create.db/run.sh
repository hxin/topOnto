#!/bin/sh
BASEDIR=$(dirname $0)

if [ -f $BASEDIR/config.txt ] ; then
. $BASEDIR/config.txt
else
	echo "config file $BASEDIR/config.txt is missing, aborting..."
	exit 1
fi

mkdir $BASEDIR/$PACKAGE_NAME
tar -zxvf topOnto.skeleton.db.tar.gz -C $BASEDIR/$PACKAGE_NAME

cat << CONFIG >$BASEDIR/$PACKAGE_NAME/inst/scripts/config
DBSCHEMA=$DBSCHEMA
DBSCHEMAVERSION=$DBSCHEMAVERSION
SOURCENAME=$SOURCENAME
SOURCURL=$SOURCURL
SOURCEDATE=$SOURCEDATE
ONLY_IS_A=$ONLY_IS_A
VERSION=$VERSION
CONFIG

cat << DESCRIPTION >$BASEDIR/$PACKAGE_NAME/DESCRIPTION
Package: $PACKAGE_NAME
Type: Package
Title: TopOnto ontology db package
Version: 0.99.0
Date: 2015-02-05
Author: Xin He
Maintainer:Xin He <xin.he@ed.ac.uk>
Depends: R (>= 2.7.0), methods, AnnotationDbi (>= 1.26.0), graph (>= 1.42.0)
Suggests: knitr
Description: TopOnto ontology db package. This package provides objects that represent an ontology.
biocViews: AnnotationData, FunctionalAnnotation
License: LGPL
VignetteBuilder: knitr
DESCRIPTION

cat << ONTBASE >$BASEDIR/$PACKAGE_NAME/man/ONTBASE.Rd
\name{$PACKAGE_NAME}
\alias{$PACKAGE_NAME}
\alias{ONT}
\title{Bioconductor annotation data package}
\description{

  Welcome to the $PACKAGE_NAME annotation Package.  The purpose of this package
  is to provide detailed information for the use to the package topOnto. 
  This package is updated biannually.
  
  You can learn what objects this package supports with the following command:
  
  \code{ls("package:$PACKAGE_NAME")} 
  
  Each of these objects has their own manual page detailing where
  relevant data was obtained along with some examples of how to use it.
}
\examples{
  ls("package:$PACKAGE_NAME")
}
\keyword{datasets}
ONTBASE

sed -i.bak 's/xxx/HDO/g' $BASEDIR/$PACKAGE_NAME/vignettes/my-vignette.Rmd | rm -rf $BASEDIR/$PACKAGE_NAME/vignettes/*.bak

cp $OBO $BASEDIR/$PACKAGE_NAME/inst/scripts/ont.obo

#echo "TO tun the script to create topOnto.db, you need to install RSQLite. DO you want to install RSQLite now?[y/n]"
#read ans

#[ -z $ans ] && ans='y'

#[ $ans = 'y' ] && (cd $BASEDIR/$PACKAGE_NAME/inst/extdata/ && R CMD INSTALL RSQLite_0.11.4.tar.gz)

(cd $BASEDIR/$PACKAGE_NAME/inst/scripts/ && chmod +x ./batch_run.sh && sh ./batch_run.sh && mv ./DB.sqlite ./../extdata/)

R CMD INSTALL --no-multiarch --with-keep.source $BASEDIR/$PACKAGE_NAME


