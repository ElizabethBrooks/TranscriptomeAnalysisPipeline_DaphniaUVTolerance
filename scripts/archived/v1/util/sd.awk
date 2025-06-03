#/usr/bin/env awk
{
	sum+=$0;
	a[NR]=$0
}
END{
	for(i in a)y+=(a[i]-(sum/NR))^2;
	print sqrt(y/(NR-1))
}