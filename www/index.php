
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

The web site <a href="http://www.stat.wisc.edu/~yandell/pda">http://www.stat.wisc.edu/~yandell/pda</a> contains free information from
<A HREF="http://www.stat.wisc.edu/~yandell/">Brian Yandell</A>'s (1997)
book <A HREF="summary.html">Practical Data Analysis for Designed
Experiments</A>.  PDA was published in January 1997 by
<A HREF="http://www.crcpress.com/www/chaphall.htm">Chapman &amp; Hall/CRC Press</A>,
<A HREF="http://www.crcpress.com/catalog/C6341.htm">ISBN 0-412-06341-7</A>. 
<ul>
<li> <A HREF="http://www.stat.wisc.edu/~yandell/pda/summary.html">Book Summary</A>
<li>	<A HREF="http://www.stat.wisc.edu/~yandell/pda/outline/">Book Outline</A>
<li> 	<A HREF="http://www.stat.wisc.edu/~yandell/pda/errata.pdf">Book Errata</A>
<li> <a href="http://www.stat.wisc.edu/~yandell/pda/data">Data Used in book (includes SAS and Splus code)</a>
<li> Publisher: <A
HREF="http://www.crcpress.com/www/chaphall.htm">Chapman &amp; Hall/CRC Press</A> (search: yandell)
     <LI> Book Reviews:
<UL>
  <LI> <a href="http://links.jstor.org/sici?sici=0040-1706%28199805%2940%3A2%3C154%3APDAFDE%3E2.0.CO%3B2-L">Grego JM (1998)
     <CITE>Technometrics 40</CITE>, 154--155.</a>
<LI> <a href="http://links.jstor.org/sici?sici=0006-341X%28199812%2954%3A4%3C1678%3APDAFDE%3E2.0.CO%3B2-K">Talbot M (1998)
     <CITE>Biometrics 54</CITE>, 1678.</a>
</UL>
</UL>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
