
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
	<title>GUTS Project<?php //echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
	<link href="./styles.css" rel="stylesheet" type="text/css" />
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
/*
if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')) {
  $contents = '';
  while (!feof($handle)) {
    $contents .= fread($handle, 8192);
  }
  fclose($handle);
  echo $contents;
}
*/
?>

<!-- end of project description -->

<h1>GUTS Project</h1>


<h2>News</h2>

<ul>
	<li>GUTS version 1.2.4 updates package documentation and improves test access. Package functionality was not changed.</li>
	<li>GUTS version 1.2.3 fixes an incorrect IT damage calculation in case of complex profiles. Since version 1.2 in few cases of specific exposure dynamics calculation of IT damage and survival estimates have gone wrong.</li>
	<li>GUTS version 1.2.2 fixes C++ code that did not compile on one computer platform supported by CRAN and improves checks of model parameter values.</li>
	<li>GUTS version 1.2.1 fixes C++ code that caused compiler warnings on some computer platforms supported by CRAN and improves consistency with C++-11 style</li>
	<li>GUTS version 1.2 introduces faster implementations for special cases GUTS-RED-IT and GUTS-RED-SD.<br>The now analytical solution for GUTS-RED-IT can be approximated by the previous version if using large M and N. The new version does not use parameters M and N any longer.<br>GUTS-RED-SD now directly implements the delta distribution and does no longer use N.<br>Further, features are: modular c++-code, template adaptation to Rcpp data structures, adaptation to ERA as suggested in EFSA PPR (2018)</li>
	<li>GUTS version 1.1.1 adds access to the sum of squares and the survival-probability prediction error.</li>
	<li>GUTS version 1.1.0 was extended to provide different distribution for individual tolerance (log-logistic, log-normal and a flexible option for applying any possible distribution function from R).</li>
	<li>New vignettes are provided documenting the process of model calibration.</li>
	<li><strong>Important:</strong> GUTS version &gt; 1.0.0 is a complete re-implementation of the former GUTS package!  Please consult the man page for more information.</li>
	<!--li>Extensive documentary material is under preparation. Updates will follow soon.</li-->
	<!-- li>The current documentation (see below) is outdated. An updated version will follow soon.</li -->
	<!--li>Current <a href="./files/GUTS_0.9.5.tgz">package as OS X binary</a>. (This may work or not.)</li-->
</ul>


<h2>Publications</h2>

<ul>
	<li>Albert, C., Vogel, S., and Ashauer, R. (2016). Computationally efficient implementation of a novel algorithm for the General Unified Threshold Model of Survival (GUTS). PLOS Computational Biology, 12(6), e1004978. doi: <a href="http://dx.doi.org/10.1371/journal.pcbi.1004978" target="_blank">10.1371/journal.pcbi.1004978</a>.</li>
	<li>EFSA PPR Panel (EFSA Panel on Plant Protection Products and their Residues), Ockleford, C., Adriaanse, P., Berny, P., Brock, T., Duquesne, S., Grilli, S., Hernandez-Jerez, A.F., Bennekou, S.H., Klein, M., Kuhl, T., Laskowski, R., Machera, K., Pelkonen, O., Pieper, S., Smith, R.H., Stemmer, M., Sundh, I., Tiktak, A., Topping, C.J., Wolterink, G., Cedergreen, N., Charles, S., Focks, A., Reed, M., Arena, M., Ippolito, A., Byers, H. and Teodorovic, I. (2018). Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms. EFSA Journal, 16(8):5377, 188 pp. doi: <a href="https://doi.org/10.2903/j.efsa.2018.5377" target="_blank">10.2903/j.efsa.2018.5377</a>.</li>
	<li>Jager, T., Albert, C., Preuss T., and Ashauer R. (2011). General Unified Threshold Model of Survival – a toxicokinetic toxicodynamic framework for ecotoxicology. Environmental Science &amp; Technology, 45(7), 2529–2540, doi: <a href="http://dx.doi.org/10.1021/es103092a" target="_blank">10.1021/es103092a</a>.</li>
</ul>


<h2>Links</h2>

<ul>
	<li>GUTS R-Forge mailing list (news and updates): <a href="https://lists.r-forge.r-project.org/mailman/listinfo/guts-users">https://lists.r-forge.r-project.org/mailman/listinfo/guts-users</a></li>

	<li><a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">Project summary page on R-Forge</a></li>

	<li>GUTS on CRAN: <a href="http://cran.r-project.org/web/packages/GUTS/index.html">http://cran.r-project.org/web/packages/GUTS/index.html</a></li>

	<li>Background of GUTS: <a href="http://www.debtox.info/about_guts.php">http://www.debtox.info/about_guts.php</a></li>

	<li>GUTS modelling: <a href="https://www.ecotoxmodels.org/guts/">https://www.ecotoxmodels.org/guts/</a></li>

	<li>EasyGUTS a user-friendly and freely available software for TK/TD modelling of survival (using the GUTS R package): <a href="https://rifcon.de/wp-content/uploads/2018/08/setac2018_nickisch_EasyGUTS.pdf">https://rifcon.de/wp-content/uploads/2018/08/setac2018_nickisch_EasyGUTS.pdf</a></li>

</ul>

<p>Last updated: 2022-02-01</p>

</body>
</html>
