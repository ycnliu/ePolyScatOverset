#
cat <<eoi >./manualdroppages/manual.html
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN" "http://www.w3.org/TR/html4/frameset.dtd">
<html>
<head>
<title>ePolyScat User's Manual</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
</head>

<frameset cols="20%,*" framespacing="1" frameborder="yes" border="1" bordercolor="#333333">
  <frame src="contents" name="contents" scrolling="NO" noresize>
  <frame src="front" name="material">
</frameset>
<noframes><body>

</body></noframes>
</html>
eoi
#
cp ./manualdroppages/manual.html ./manualdroppages/0manual.html
#
cat <<eoi >./manualdroppages/front.html
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">

<HTML>
<HEAD>
<TITLE>First Page</TITLE>
</HEAD>

<BODY LINK=navy VLINK=navy BGCOLOR="#FFFFCC">
<BR>

<!--Table of Contents-->
<center>
<br />
<font size="7" color="purple">ePolyScat</font><br />
<h2>Date: ${GITDATE}</h2>
<h2>Commit: ${GITSHA1:0:7}</h2>
<br /><br />
<h2>R. R. Lucchese, N. Sanna, A. P. P. Natalense, and F. A. Gianturco </h2>
<br />
<p>
This is a series of programs for computing electron-molecule scattering
and molecular photoionization cross sections.
</p>
<br />
<br />
<br />
<br />
<p>
Send questions, comments, and requests for the programs to
<a href=mailto:rlucchese@lbl.gov> rlucchese@lbl.gov </a>
</p>

<p>
These programs are installed and can be used on the 
<a href="https://amosgateway.org" target="_blank"> AMOS gateway.
</p>
</center>

</BODY>
</HTML>
eoi
#
cat <<eox >./manualdroppages/DiffE3to2020.html
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
<head>
<title>DiffE3to2020</title>
</head>

<body link="navy" vlink="navy" bgcolor="FFFFCC">

<h2>Added commands and data records between the E3 version and 2020 </h2>
<p>
<code>AbSym</code> - do calculation with largest abelian subgroup
</p>
<p>
<code>AtomSym</code> - pick the symmetry to use when calculations with an atom
</p>
<p>
<code>Convert</code> <var>FileName</var> <code>'GAMESS'</code> - input filter for the US GAMESS quantum chemistry program
</p>
<p>
<code>DumpMesa</code> - write basis set and orbital information for the MEAS codes
</p>
<p>
<code>EngToler</code> - used in vibrational averaging programs
</p>
<p>
<code>FegeScale</code> - factor to scale the approximate local exchange potential
</p>
<p>
<code>GenFormScat</code> - generate scattering potential formulas for electron scattering
</p>
<p>
<code>MFDCS</code> - molecular frame differential scattering for electron-molecule scattering include scattering from a charged target
</p>
<p>
<code>NoSym</code> - do calculation with C1 symmetry
</p>
<p>
<code>PlotDataGrid</code> - controls output to the FileName 'PlotData' file
</p>
<p>
<code>PrintBlm</code> - print symmetry adapted angular functions
</p>
<p>
<code>RotOrientAsym</code> - compute MFPAD including rotation between ionization and fragmentation
</p>
<p>
<code>SchmidtOrth</code> - Schmidt orthogonalization of orbitals that have been expanded in ExpOrb
</p>
<p>
<code>ViewOrbPartialWave</code> - define a partial wave grid to use with ViewOrb
</p>
</html>
eox
$bindir/MakeManualDropPages.exe <<eoi
 './manualdroppages'
 '.' 'Intro' 'README2' ' ' 'Introduction'
 './src' 'Source', 'ePolyScat.f90' '! **' 'Commands'
 './src' 'Source' 'GetDataRecordDef.f90' '! **' 'Data Records'
 'created' 'SymmetryLabels.html' 'Symmetry' ' ' 'Symmetry Labels'
 './tests' 'Tests' 'README' ' ' 'Sample Jobs'
 '.' 'Makefile' 'README3' ' ' 'Makefile Options'
 './src' 'Source' 'README2' ' ' 'Environment Variables'
 'created' 'UtilityProgs.html' 'Utiltiy' ' ' 'Utility Programs'
eoi
#
#  './src' 'Source' 'README' ' ' 'Program Inputs'
chmod a+r manualdroppages/*.html
#
