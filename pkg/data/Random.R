genRandom = function(sd=c(20,5),size=c(10,10),mean=100)
{
  Random = data.frame( matrix( rnorm( prod( size ), 0, sd[2] ),
                               size[2], size[1] ))
  names(Random) = seq( size[1] )
  for (i in names(Random))
    Random[[i]] = Random[[i]] + rnorm( 1, mean, sd[1] )
  Random
}
Random = genRandom()
rm( genRandom )
