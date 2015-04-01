
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php
if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')) {
  $contents = '';
  while (!feof($handle)) {
    $contents .= fread($handle, 8192);
  }
  fclose($handle);
  echo $contents;
}
?>

<!-- end of project description -->

<p>The <strong>project summary page</strong> can be found <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>.
</p>

<p>A <strong>documentation</strong> can be found <a href="http://guts.r-forge.r-project.org/doc/GUTS Documentation.pdf"><strong>here</strong></a>.
</p>

<p>For questions and news on updates you can subscribe to the mailing list here: https://lists.r-forge.r-project.org/mailman/listinfo/guts-users
</p>

<p>To learn more about the background of GUTS, please refer to one of the following web sites:
<ul>
  <li><a href="http://www.debtox.info/about_guts.php">http://www.debtox.info/about_guts.php</a></li>
  <li><a href="http://www.ecotoxmodels.org/research-publications-projects/toxicokinetic-toxicodynamic-models/">http://www.ecotoxmodels.org/research-publications-projects/toxicokinetic-toxicodynamic-models/</a></li>
  <li><a href="http://dx.doi.org/10.1021/es103092a">http://dx.doi.org/10.1021/es103092a</a></li>
</ul>

<p>Last updated: 2015-04-01</p>

</body>
</html>
