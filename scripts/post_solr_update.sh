# post.jar is included in solr package

BASEDIR=$(dirname $0)
POSTJAR=$BASEDIR/post.jar
URL=$KBASE_SOLR/$1/update

java -Durl=$URL -Dauto=yes -Dfiletypes=json -Dcommit=yes -Dout=yes -Drecursive=yes -jar $POSTJAR $2
