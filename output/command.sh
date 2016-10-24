#!/bin/bash
echo '#!/bin/bash'
echo 'rm pendulum.gif 2> /dev/null'
for i in `seq 0 499`;
do
	echo 'convert pgplot'$i'.png -resize 680x680\! pgplot'$i'.png'
done
echo -n 'convert -loop 0 '
for i in `seq 0 499`;
do
	echo -n pgplot$i.png' '
done
echo 'pendulum.gif'