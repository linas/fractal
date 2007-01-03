#! /usr/bin/env perl

for ($i=0; $i<10; $i++)
{
	# system ("ls ptest-$i-*.flo ");
	system ("cat phase-small-$i-*.flo | /home/linas/src/fractal/image/flo2mtv |mtvtoppm | pnmtopng > msmall-$i.png");
}

