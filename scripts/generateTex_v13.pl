#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

sub iso8601_date {
  die unless $_[0] =~ m/^(\d\d\d\d)-(\d\d)-(\d\d)$/;
  return  ($1*365 + $2*30+ $3);
}

#Options Variablen setzen
my $in_var_fh         = "summary.csv";
my $in_mut_fh         = "summary.csv";
my $syn_fh            = "synopsis.tex";
my $mutations_fh      = "VaQuERo/resources/mutations_list_grouped_pango_codonPhased_2024-02-22_Europe_adapted.csv";
my $smutations_fh     = "VaQuERo/resources/mutations_special_2022-12-21.csv";
my $groups_fh         = "VaQuERo/resources/groupMembers_pango_codonPhased_2023-11-11_Europe.csv";
my $fig_fh            = "Mutations_of_Interest.pdf";
my @fig               = ();
my $fig               = "";
my @fig_head          = ();
my $fig_head          = "";
my $tab_fh            = "growthrate_per_mutation_pois.csv";
my $tab               = "RBD Mutationen mit Wachstum in mehr als 2 Kläranlagen.";
my $resulttable_fh    = "output-general/globalFittedData.csv";


my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $date = "${year}-${mon}-${mday}";

GetOptions(
  "iv:s"      => \$in_var_fh,
  "im:s"      => \$in_mut_fh,
  "s:s"       => \$syn_fh,
  "d:s"       => \$date,
  "m:s"       => \$mutations_fh,
  "n:s"       => \$smutations_fh,
  "t:s"       => \$resulttable_fh,
  "g:s"       => \$groups_fh,
  "f:s"       => \$fig_fh,
  "h=s{1,}"   => \@fig_head,
  "v=s{1,}"   => \@fig,
  "ta:s"      => \$tab_fh,
  "tat:s"     => \$tab
);

print STDERR "LOG: input variants  == $in_var_fh \n";
print STDERR "LOG: input mutation  == $in_mut_fh \n";
print STDERR "LOG: date   == $date \n";

my %data = ();
my %data_m = ();
my %mutation = ();
my %smutation = ();
my %restable     = ();
my %groups     = ();

print STDERR "$date\n";
my $refdate_decimal = iso8601_date($date);

open T, "< $resulttable_fh" or die "Error: can t open $resulttable_fh: $!\n";
while(<T>){
  chomp;
  my @F = split"\t", $_;
  next if ($F[0] eq "variant");
  next unless ($F[4]);
  if($F[4]){
    my $cmpdate_decimal = iso8601_date($F[4]);
    unless(($refdate_decimal-$cmpdate_decimal) < 31){
      next;
    }
  }

  if($F[5] > 0){
    $restable{$F[1]}->{$F[2]}->{$F[4]}->{$F[0]} = $F[5];
  }
}
close T;

open M, "< $mutations_fh" or die "Error: can t open $mutations_fh: $!\n";
while(<M>){
  chomp;
  my @F = split",", $_;
  next if ($F[0] eq "Variants" && $F[1] eq "Chromosom");
  my @v = split";", $F[0];
  my @s = split";", $F[6];
  my $aa = "";
  if ($F[7]){
    my @aas = split(":", $F[7]);
    $aa = $aas[-1];
  }
  my $ref = $F[3];
  my $alt = $F[4];
  my $nuc = $F[8];
  if (length($ref) > 4 || length($alt) > 4){
    $nuc = "longIndel";
    if (length($ref) > length($alt)){
      $ref = "del:".abs(length($ref)-length($alt));
    }
    elsif(length($ref) < length($alt)){
      $alt = "ins:".abs(length($ref)-length($alt));
    }
  }
  my $count = scalar(@v);
  foreach my $i (0..$#v){
    my $v = $v[$i];
    my $s = $s[$i];
    $mutation{$v}->{$F[1]}->{$F[2]}->{$ref}->{$alt}->{$F[5]}->{$s}->{$aa}->{$nuc} = $count;
  }
}
close M;


open SM, "< $smutations_fh" or die "Error: can t open $smutations_fh: $!\n";
while(<SM>){
  chomp;
  my @F = split",", $_;
  next if ($F[0] eq "Variants" && $F[1] eq "Chromosom");
  my @v = split";", $F[0];
  my @s = split";", $F[6];
  my $aa = "";
  if ($F[7]){
    my @aas = split(":", $F[7]);
    $aa = $aas[-1];
  }
  my $ref = $F[3];
  my $alt = $F[4];
  my $nuc = $F[8];

  if (length($ref) > 4 || length($alt) > 4){
    $nuc = "longIndel";
    if (length($ref) > length($alt)){
      $ref = "del:".abs(length($ref)-length($alt));
    }
    elsif(length($ref) < length($alt)){
      $alt = "ins:".abs(length($ref)-length($alt));
    }
  }
  my $count = scalar(@v);
  foreach my $i (0..$#v){
    my $v = $v[$i];
    my $s = $s[$i];
    $smutation{$v}->{$F[1]}->{$F[2]}->{$ref}->{$alt}->{$F[5]}->{$s}->{$aa}->{$nuc} = $count;
  }
}
close SM;

open G, "< $groups_fh" or die "Error: can t open $groups_fh: $!\n";
while(<G>){
  chomp;
  my @F = split"\t", $_;
  next if ($F[1] eq "variants" && $F[2] eq "groupName");

  $groups{$F[2]}->{$F[1]}++;
}
close G;


open I, "< $in_var_fh" or die "Error: can t open plot2path (variants) file $in_var_fh: $!\n";
while(<I>){
  chomp;
  my @F = split"\t", $_;
  if(scalar(@F) <= 1){
    print STDERR "Warning: less than 2 entries in line $.\n";
  }
  elsif(scalar(@F) == 2){
    $data{$F[0]} = $F[1];
  }
  elsif(scalar(@F) == 3){
    $data{$F[0]}->{$F[1]} = $F[2];
  }
  elsif(scalar(@F) == 4){
    $data{$F[0]}->{$F[1]}->{$F[2]} = $F[3];
  }
  elsif(scalar(@F) == 5){
    $data{$F[0]}->{$F[1]}->{$F[2]}->{$F[3]} = $F[4];
  }
  else{
    print STDERR "Warning: more than 5 entries in line $.\n";
  }
}
close I;

open I, "< $in_mut_fh" or die "Error: can t open plot2path (mutations) file $in_mut_fh: $!\n";
while(<I>){
  chomp;
  my @F = split"\t", $_;
  if (scalar(@F) != 4){
    next;
  }

  if(scalar(@F) == 4 && $F[0] eq "cuadrilla") {
    $data_m{$F[0]}->{$F[1]} = $F[3];
  } elsif(scalar(@F) == 4 && $F[0] eq "excessmutations"){
    $data_m{$F[0]}->{$F[1]} = $F[3];
  }
}
close I;


# print header + external Synopsis text
&printHeader($syn_fh);

# print subsections for maps
if (defined($data{statusmapplot})){
  my $path = $data{statusmapplot}->{status};
  &printStatusMap($path);
}
else{
  print STDERR "Warning: no overview map\n";
}

# print subsections for maps
if (defined($data{mapplot})){
  foreach my $variant ( sort keys %{$data{mapplot}} ){
    my $path = $data{mapplot}->{$variant};
    &printMapPerVariant($variant, $path);
  }
}
else{
  print STDERR "Warning: no overview map\n";
}

if(-f "$fig_fh"){
  $fig = join " ",@fig;
  $fig_head = join " ", @fig_head;

  &printFig($fig_fh, $fig, $fig_head);
  print STDERR "Log: use plot <$fig_fh> as special Figure\n";
  print STDERR "Log: use section heading <$fig_head> for special Figure\n";
  print STDERR "Log: use caption <$fig> for special Figure\n";
}
else{
  print STDERR "Warning: no plot found by name $fig_fh\n";
}

# print subsections for growing excess mutation overview
if (defined($data_m{excessmutations}->{kineticPlot}) && defined($data_m{excessmutations}->{outbreakInfoPlot}) && defined($data_m{excessmutations}->{table}) && defined($data_m{cuadrilla}->{plot})) {
    my $path_plot = $data_m{excessmutations}->{kineticPlot};
    my $path_outbreakinfo = $data_m{excessmutations}->{outbreakInfoPlot};
    my $path_table = $data_m{excessmutations}->{table};
    my $path_plot_cuadrilla = $data_m{cuadrilla}->{plot};
    &printMutationOverview($path_plot, $path_outbreakinfo, $path_table, $path_plot_cuadrilla);
}
else{
  print STDERR "Warning: no growing_excessmutations plots\n";
}


#print header section kläranlagen ergebnisse
print '\clearpage'."\n";
print '\section{Ergebnisse per Kläranlage}'."\n";
print '\subsection{Ergebnisse per Kläranlage mit Daten aus aktuellem Run}'."\n";
if (defined($data{variantDetail}->{current}) && defined($data{stackOverview}->{current})){
  foreach my $wwplant ( sort keys %{$data{variantDetail}->{current}} ){
    my $wwplant_tidy = $wwplant;
    $wwplant_tidy =~ s/_/ /g;
    print '\clearpage'."\n";
    print '\subsubsection{'."$wwplant_tidy".'}'."\n";

    # print graphical synapsis
    if ( defined($data{sankey}->{WWTP}->{$wwplant}) ){
      my $path = $data{sankey}->{WWTP}->{$wwplant};
      &printSankey($path, $wwplant_tidy);
    }

    if( defined($data{stackOverview}->{current}->{$wwplant}->{all}) ){
      my $path1 = $data{stackOverview}->{current}->{$wwplant}->{all};
      &printStackPlot($wwplant_tidy, $path1);
    }

    if ( defined($data{variantDetail}->{current}->{$wwplant}->{all}) ){
      my $path2 = $data{variantDetail}->{current}->{$wwplant}->{all};
      #&printDetailPlot($wwplant_tidy, $path2);
    }
  }
}


if(0){
  print '\clearpage'."\n";
  print '\subsection{Ergebnisse per Kläranlage ohne Daten aus aktuellem Run}'."\n";
  if (defined($data{variantDetail}->{old}) && defined($data{stackOverview}->{old})){
    foreach my $wwplant ( sort keys %{$data{variantDetail}->{old}} ){
      my $wwplant_tidy = $wwplant;
      $wwplant_tidy =~ s/_/ /g;
      print '\clearpage'."\n";
      print '\subsubsection{'."$wwplant_tidy".'}'."\n";

      if( defined($data{stackOverview}->{old}->{$wwplant}->{all}) ){
        my $path1 = $data{stackOverview}->{old}->{$wwplant}->{all};

        &printStackPlot($wwplant_tidy, $path1);
      }
      if ( defined($data{variantDetail}->{old}->{$wwplant}->{VoI}) ){
        my $path2 = $data{variantDetail}->{old}->{$wwplant}->{VoI};

        #&printDetailPlot($wwplant_tidy, $path2);
      }
      if ( defined($data{specialMutations}->{old}->{$wwplant}->{VoI}) ){
        my $path3 = $data{specialMutations}->{old}->{$wwplant}->{VoI};

        # &printSpecialMutPlot($wwplant_tidy, $path3);
      }
    }
  }
}

# print results of last 30 days as table
&printResults();

# print Variant Group Appendix
&groupsAppendix();

# print Marker Mutations Appendix
&mutationsAppendix();

# print Special Mutations Appendix
if($smutations_fh){
&smutationsAppendix();
}

# print Daten und Methodik und Haftungsausschuss
&printDMH();

# print footer
&printFooter();


#print Dumper(%data);
####

sub printGraphicalSynopsis{
my $fh = shift;
my $lbl = shift;
my $txt='

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$fh".'}
    \caption{Graphische Darstellung der Varianten Quantifizierung für die neusten Proben jeder Kläranlage '."$lbl".' für den Berichtszeitraum. Dargestellte Werte entsprechende Mittelwerten, gewichtet nach Bevölkerungsgröße. Phylogenetische Verwandtschaft wird durch die Baumorganisation visualisiert.}
  \end{center}
\end{figure}
';

print $txt;
}

sub printGraphicalSynopsis2{
my $fh = shift;
my $lbl = shift;
my $txt='

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$fh".'}
    \caption{'."$lbl".'}
  \end{center}
\end{figure}
';

print $txt;
}


sub printSankey{
my $fh = shift;
my $lbl = shift;
my $txt='
\FloatBarrier
{\bf {\large Überblick aller detektierten Varianten in der aktuellen Probe}}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$fh".'}
    \caption{Graphische Darstellung der Varianten Quantifizierung für die neusten Proben der Kläranlage '."$lbl".'. Phylogenetische Verwandtschaft wird durch die Baumorganisation visualisiert.}
  \end{center}
\end{figure}
';

print $txt;
}


sub printFig{
my $fh = shift;
my $lbl = shift;
my $head = shift;
my $txt='
\clearpage
\section{Fokus Analyse}
\subsection{'."$head".'}
\label{Fokus}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=.66\textwidth]{'."$fh".'}
    \caption{'."$lbl".'}
  \end{center}
\end{figure}
';

print $txt;
}

sub printTab{

my $fh = shift;
my $headline = shift;


my $txt = '

\subsection{Einzel Mutationen Kinetik Tabular}

'."$headline".'

\begin{footnotesize}
\begin{longtable}{p{.22\textwidth} | p{.15\textwidth} | p{.15\textwidth} | p{.15\textwidth} | p{.25\textwidth}}
\hline
Mutation & Mittleres R $w^{-1}$ & Max R $w^{-1}$ & Nr. WWTP & cov-spectrum \\\\
\hline
';

open D, "< $fh" or die "Can't open $fh: $!\n";
while(<D>){
  chomp;
  my @F = split"\t", $_;
  next if ($F[0] eq "mutation");
  my $label = $F[1];
  my $median =sprintf("%.3f", $F[3]);
  my $max =sprintf("%.3f", $F[4]);
  my $n = $F[7];
  next if ($n <= 2);
  my $link = join("", '\\href{https://cov-spectrum.org/explore/Austria/AllSamples/Past2M/variants?nextcladeCoverageFrom=0.5&nextcladeCoverageTo=1&nucMutations=', $F[0],'&}{', $F[0], '}');
  $txt .= "
$label & $median & $max & $n & $link \\\\
   ";
}

$txt .= '
\label{tab:rbd}
\end{longtable}
\end{footnotesize}
';

$txt =~ s/_/\\_/g;
print $txt;
}

sub printFooter {
my $txt = '


%----------------------------------------------------------------------------------------
%	BIBLIOGRAPHY
%----------------------------------------------------------------------------------------

\bibliographystyle{apalike}

\bibliography{sample}

%----------------------------------------------------------------------------------------


\end{document}
';
print $txt;
};

sub printResults{

my $txt = '
\clearpage

\section{Resultate in Tabellen Form}

Tabelle mit allen abgeleiteten Varianten Häufigkeiten für die letzten 21 Tage.

\begin{footnotesize}
\begin{longtable}{p{.20\textwidth} | p{.33\textwidth} | p{.10\textwidth} | p{.10\textwidth} | p{.08\textwidth} }
\hline
Location ID & Location Name & Datum & Variante & Anteil\\\\
\hline
';

foreach my $ln (sort keys %restable){
  foreach my $li (sort keys %{$restable{$ln}}){
    foreach my $d (sort keys %{$restable{$ln}->{$li}}){
      foreach my $v (sort keys %{$restable{$ln}->{$li}->{$d}}){
         my $value =sprintf("%.1f \\%%", 100*$restable{$ln}->{$li}->{$d}->{$v});
         $txt .= "
$ln & $li & $d & $v & $value \\\\
          ";
      }
      $txt .= '\hdashline
      ';
    }
  }
  $txt .= '\hline \\\\
  \hline
      ';
}
$txt .= '
\label{tab:results}
\end{longtable}
\end{footnotesize}
';

$txt =~ s/_/\\_/g;
print $txt;
}

sub mutationsAppendix{

my $txt = '
\clearpage
\section{Appendix: Verwendete Markermutations}

Auflistung aller verwendeter Markermutationen für die untersuchten Varianten. Varianten die hier nicht aufgeführt sind können in den Abwasserproben auch nicht quantifiziert werden. Die Auswahl der Varianten berücksichtigt zum einen, deren Einstufung von ECDC, die Häufigkeit in Österreich laut GISAID, und deren globalen Wachstum in den letzten 6 Monaten laut cov-spectrum. Die Spalte Eindeutigkeit gibt an für wieviele betrachteten Varianten die entsprechende Mutation als sensitive Markermutation definiert wurde. Spezifische Markermutationen ("unique marker") sind mit einer 1 in dieser Spalte markiert. Sie Spalte Sensitivität gibt den Anteil der Samples an, die in Pango-Designation die dieser Variante zugeordnet werden, und diese Mutation tragen. Lange INDELS wurden in der Tabelle, zugunsten der Lesbarkeit, gekürzt angegeben.

\begin{footnotesize}
\begin{longtable}{p{.08\textwidth} | p{.05\textwidth} | p{.06\textwidth} | p{.06\textwidth} | p{.10\textwidth} | p{.10\textwidth} | p{.11\textwidth} | p{.11\textwidth} | p{.03\textwidth}}
\hline
Variante & Position & REF & ALT & Gen & Sensitivität & AA & NUC & Eindeutigkeit \\\\
\hline
';

foreach my $var (sort keys %mutation){
  foreach my $chrom (sort keys %{$mutation{$var}}){
    foreach my $pos (sort keys %{$mutation{$var}->{$chrom}}){
      foreach my $ref (sort keys %{$mutation{$var}->{$chrom}->{$pos}}){
        foreach my $alt (sort keys %{$mutation{$var}->{$chrom}->{$pos}->{$ref}}){
          foreach my $gen (sort keys %{$mutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}}){
            foreach my $sens (sort keys %{$mutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}}){
              foreach my $aa (sort keys %{$mutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}->{$sens}}){
                foreach my $nuc (sort keys %{$mutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}->{$sens}->{$aa}}){
                  my $uniq = $mutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}->{$sens}->{$aa}->{$nuc};

                  $txt .= "
$var & $pos & $ref & $alt & $gen & $sens & $aa & $nuc & $uniq \\\\
                  ";
                }
              }
            }
          }
        }
      }
    }
  }
}
$txt .= '
\label{tab:markermuts}
\end{longtable}
\end{footnotesize}
';

$txt =~ s/_/\\_/g;
print $txt;
}


sub smutationsAppendix{

my $txt = '
\clearpage
\section{Appendix: Berücksichtigte Spezial-Mutationen}
\label{a:specmut}

Für manche Pangolin Varianten reichen die spezifischen Markermutationen nicht um sie mit dem hier verwendeten in VaQuERo implementierten Ansatz zu unterscheiden. Um trotzdem einen Eindruck gewinnen zu können werden selektive zusätzlich Spezial-Mutationen für gewisse Varianten mit ihrer Mutationshäufigkeit graphisch dargestellt. Folgende Tabelle listet alle verwendeten Spezial-Mutationen und deren assoziierte Variante.

\begin{footnotesize}
\begin{longtable}{p{.16\textwidth} | p{.08\textwidth} | p{.06\textwidth} | p{.06\textwidth} | p{.10\textwidth} | p{.11\textwidth} | p{.11\textwidth} }
\hline
Ass. Variante & Position & REF & ALT & Gen & AA & NUC \\\\
\hline
';

foreach my $var (sort keys %smutation){
  foreach my $chrom (sort keys %{$smutation{$var}}){
    foreach my $pos (sort keys %{$smutation{$var}->{$chrom}}){
      foreach my $ref (sort keys %{$smutation{$var}->{$chrom}->{$pos}}){
        foreach my $alt (sort keys %{$smutation{$var}->{$chrom}->{$pos}->{$ref}}){
          foreach my $gen (sort keys %{$smutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}}){
            foreach my $sens (sort keys %{$smutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}}){
              foreach my $aa (sort keys %{$smutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}->{$sens}}){
                foreach my $nuc (sort keys %{$smutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}->{$sens}->{$aa}}){
                  my $uniq = $smutation{$var}->{$chrom}->{$pos}->{$ref}->{$alt}->{$gen}->{$sens}->{$aa}->{$nuc};

                  $txt .= "
$var & $pos & $ref & $alt & $gen & $aa & $nuc \\\\
                  ";
                }
              }
            }
          }
        }
      }
    }
  }
}
$txt .= '
\label{tab:smarkermuts}
\end{longtable}
\end{footnotesize}
';
$txt =~ s/_/\\_/g;
print $txt;
}

sub groupsAppendix{

my $txt = '
\clearpage
\section{Appendix: Varianten}
\label{appendix:Varianten}

Im vorliegenden Abwasser-Varianten Bericht werden Virus Varianten nach der gängigen
Klassifizierung von \href{https://cov-lineages.org/lineage_list.html}{Pangolin}
unterschieden. Diese Klassifizierung ist sehr volatil. Manche beschriebene
Sub-Varianten unterscheiden sich in so wenigen Mutationen dass eine Unterscheidung
mit der hier angewandten Methodik nicht möglich ist. Daher werden gewisse Sub-Varianten
nicht weiter differenziert sondern aggregiert ausgegeben. Folgende Tabelle listet
alle aggregierten Varianten, falls vorhanden deren Alias und deren im Bericht
verwendete Bezeichnung.

\begin{footnotesize}
\begin{longtable}{p{.2\textwidth} | p{.2\textwidth}}
\hline
Varianten Name & Inkl. Varianten \\\\
\hline
\hline
';

foreach my $name (sort keys %groups){
  $txt .= "
\\hline
$name & $name \\\\
  ";
  foreach my $var (sort keys %{$groups{$name}}){
    next if ($name eq $var);
    $txt .= "
      & $var \\\\
  " ;
  }
}
$txt .= '
\label{tab:groups}
\end{longtable}
\end{footnotesize}
';

print $txt;
}



sub printDetailPlot{
  my $location = shift;
  my $path     = shift;

my $txt = '
\FloatBarrier
{\bf {\large Detailansicht ausgewählter Varianten}}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$path".'}
    \caption{Abgeleitete relative Häufigkeit aller detektierten VoC/VoI Variante(n) und der gemessenen Allelfrequenz der dazugehörenden uniquen Markermutationen für die gesamte Messzeitreihe für das Klärwerk '."$location".'.}
  \end{center}
\end{figure}
';
print $txt;
}

##XX##

sub prinExcessGrowPlot{
  my $location = shift;
  my $pathp     = shift;
  my $patht     = shift;

  open TT, "< $patht " or die "can t open patht ($patht): $!\n";
  my $table_content = "";
  while(<TT>){
    $table_content.=$_;
  };
  close TT;


my $txt = '
\FloatBarrier
{\bf {\large Überschuss Mutationen}}

'."$table_content".'
\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$pathp".'}
    \caption{Überschuß-Mutationen für das Klärwerk '."$location".', d.h., Mutationen die nicht durch detektierte Varianten erklärt werden können und ein wöchentliches Wachstum größer 0.02 zeigen.}
  \end{center}
\end{figure}
';
print $txt;
}

##YY##
sub printSpecialMutPlot{
  my $location = shift;
  my $path     = shift;

my $txt = '
\FloatBarrier
{\bf {\large Detailansicht ausgewählter Spezial-Mutationen}}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$path".'}
    \caption{Abgeleitete relative Häufigkeit aller detektierten VoC/VoI Variante(n) und der gemessenen Allelfrequenz der ausgewählter Spezial-Mutationen (siehe Appendix \ref{a:specmut}) für die gesamte Messzeitreihe für das Klärwerk '."$location".'.}
  \end{center}
\end{figure}
';
print $txt;
}


sub printStackPlot{
  my $location = shift;
  my $path     = shift;

my $txt = '
\FloatBarrier
{\bf {\large Überblick aller detektierten Varianten}}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.66\textwidth]{'."$path".'}
    \caption{Relative Häufigkeit aller detektierten Variante(n) für die gesamte Messzeitreihe für das Klärwerk '."$location".'.}
  \end{center}
\end{figure}
';
print $txt;
}



sub printStatusMap{
  my $path    = shift;

my $txt = '
\clearpage
\subsection{Übersicht '."Sequenzierstatus".'}
\label{map'."status".'}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.5\textwidth]{'."$path".'}
    \caption{Zusammenfassung des Sequenziererfolgs aller Proben aus der letzten Sequenziertranche. Wenn mindestens 40\% des Genomes erfolgreich rekonstruiert werden konnte, gilt die Probe als erfolgreich sequenziert und wird im Weiteren auch analysiert. Wenn weniger als 40\% aber mehr als 5\% des Genomes rekonstuiert werden konnte, wird die Probe mit "Virus detektiert" makiert, aber nicht versucht verschiedene Varianten zu quantifizieren. Darunter gilt die Probe als fehlgeschlagen. Da manche Kläranlagen mehrfach pro Woche beprobt werden, können auch fehlgeschlagene und erfolgreiche Beprobungen für die selbe Anlage in einer Sequenziertranche vorkommen. .}
  \end{center}
\end{figure}
';
print $txt;
}

sub printMutationOverview{
  my $plot  = shift;
  my $outbreakInfoPlot  = shift;
  my $table = shift;
  my $plot_cuadrilla = shift;
  open OT, "< $table" or die "can t open $table: $!\n";
  my $table_content = "";
  while(<OT>){
    $table_content.=$_;
  };
  close OT;

my $txt = '
\clearpage
\section{Übersicht wachsende Überschuss-Mutationen}
\label{excess_growth_overview}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.44\textwidth]{'."$plot".'}
    \includegraphics[width=0.44\textwidth]{'."$outbreakInfoPlot".'}
    \caption{Mutationen und deren potentielle assoziierten Varianten, die in der aktuellsten Woche österreichweit zunehmende geographische Ausbreitung und zunehmende relative Häufigkeit zeigen.}
  \end{center}
\end{figure}
'."$table_content".'

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.95\textwidth]{'."$plot_cuadrilla".'}
    \caption{Mutationen und deren potentielle assoziierten Varianten, die im betrachteten Zeitraum zumindest in zwei aufeinander folgenden Wochen österreichweit zunehmende geographische Ausbreitung und zunehmende relative Häufigkeit zeigen.}
  \end{center}
\end{figure}
';



print $txt;
}

sub printMapPerVariant{
  my $variant = shift;
  my $path    = shift;

my $txt = '
\clearpage
\subsection{Übersicht '."$variant".'}
\label{map'."$variant".'}

\begin{figure}[htpb]
  \begin{center}
    \includegraphics[width=0.5\textwidth]{'."$path".'}
    \caption{Der Farbton zeigt den relativen Anteil, falls größer Null, der Variante '."$variant".' zum Zeitpunkt der letzten Probennahme. Stellen mit einem Anteil von Null werden nicht gezeichnet.}
  \end{center}
\end{figure}

';
print $txt;
}

sub printHeader{
my $extTxtFh = shift;
my $txt1 = '
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% University/School Laboratory Report
% LaTeX Template
% Version 3.1 (25/3/14)
%
% This template has been downloaded from:
% http://www.LaTeXTemplates.com
%
% Original author:
% Linux and Unix Users Group at Virginia Tech Wiki
% (https://vtluug.org/wiki/Example_LaTeX_chem_lab_report)
%
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------------------------------
%	PACKAGES AND DOCUMENT CONFIGURATIONS
%----------------------------------------------------------------------------------------

\documentclass[a4paper,10pt]{article}

\usepackage{graphicx} % Required for the inclusion of images
\usepackage{natbib} % Required to change bibliography style to APA
\usepackage{amsmath} % Required for some math elements
\usepackage[section]{placeins}
\usepackage{titleps}
\usepackage{hyperref}
\usepackage[margin=1in]{geometry}
\usepackage{longtable}
\usepackage{arydshln}
\usepackage{caption}

\setlength\parindent{0pt} % Removes all indentation from paragraphs

\renewcommand{\labelenumi}{\alph{enumi}.} % Make numbering in the enumerate environment by letter rather than number (e.g. section 6)

\usepackage{times} % Uncomment to use the Times New Roman font

%----------------------------------------------------------------------------------------
%	DOCUMENT INFORMATION
%----------------------------------------------------------------------------------------
% header and footer definition
\newpagestyle{main}{
  \setheadrule{.3pt}% Header rule
  \sethead{} % left
          {} % center
          {\thesubsection\ \subsectiontitle} %right
  \setfootrule{.3pt}% Header rule
  \setfoot{} % left
          {} % center
          {CeMM - preliminary data} %right
%         {Vorläufige Datenauswertung: CeMM - Sars_Cov-2 Abwasser-Varianten-Monitor} % rigth
}
\pagestyle{main}

% define if sections (=1), subsection (=2) should be included in ToC
\setcounter{tocdepth}{3}

\title{Österreichischer Abwasser Sars-CoV-2 Varianten Bericht} % Title

\author{Fabian \textsc{Amman} \and Lukas \textsc{Endler} \and Anna \textsc{Schedl} \and Matthew \textsc{Thornton} \and Anastasiia \textsc{Antielieieva} \and Bethany \textsc{Dearlove} \and Andreas \textsc{Bergthaler}} % Author name

\date{\today} % Date for the report

\begin{document}

\maketitle % Insert the title, author and date
\begin{center}
  \begin{tabular}{l l }
    Kontakt: & \href{mailto:andreas.bergthaler@meduniwien.ac.at}{andreas.bergthaler@meduniwien.ac.at} \\\\
             & \href{mailto:FAmman@cemm.at}{FAmman@cemm.at} \\\\
    Berichtsdatum: & '."$date".' \\\\ % Date the experiment was performed
    Auswahlkriterium Varianten: & Berücksichtigte Varianten aus drei Quellen \\\\
                                & \\\\
                                & Alle \emph{Variants~of~concern}, \emph{of~interest}, und \emph{under~monitoring}. \\\\
                                & Definition folgt der Angabe von \href{https://www.ecdc.europa.eu/en/covid-19/variants-concern}{ECDC} und \href{https://data.who.int/dashboards/covid19/variants}{WHO} (Stand 15. Juni., 2025)\\\\
                                & \\\\
                                & Alle Varianten die laut \href{https://gisaid.org/}{GISaid} Daten in Europa in zumindest\\\\
                                & einem Land und einem Monat eine relative Häufigkeit $>$60\% aufwiesen, oder\\\\
                                & einem Land und einem Monat eine relative Häufigkeit $>$10\% seit März 2024 aufwiesen, \\\\
%                               & \\\\
%                               & Alle Varianten die aktuell laut \href{https://cov-spectrum.org/collections/24}{CoV-spectrum Collection 24}\\\\
%                               & ein Wachstumsvorteil $>$30\% (lower CI) zeigen. Detailgrad der \\\\
% 			                        & Untervarianten nach individueller Einschätzung angepasst. \\\\
                                & \\\\
                                & Siehe Varianten Liste im Appendix~\ref{appendix:Varianten} auf Seite~\pageref{appendix:Varianten} für \\\\
                                & eine detailierte Auflistung der berücksichtigten Varianten.\\\\
  \end{tabular}
\end{center}
\newpage

%----------------------------------------------------------------------------------------
%	SECTION 1
%----------------------------------------------------------------------------------------

\section{Synopsis}
';

open TXT, "< $extTxtFh" or die "can't open since $extTxtFh not readable: $!\n";
my @txt2 =<TXT>;
my $txt2 = join("", @txt2);

my $txt3 = '


\tableofcontents
\clearpage

%----------------------------------------------------------------------------------------
%	SECTION 2
%----------------------------------------------------------------------------------------

\section{Übersicht}

';
print $txt1;
print $txt2;
# print graphical synapsis
if (defined($data{sankey}->{Overview}->{Austria}) || defined($data{sankey}->{allStates}->{Austria})){
    print '\subsection{Graphische Synopsis}'."\n";
    print '\label{graphsyn}'."\n";

    if (defined($data{stacked}->{Overview}->{Austria})){
      my $path = $data{stacked}->{Overview}->{Austria};
      my $text = "Relative Varianten Häufigkeit über den gesamten Pandemieverlauf, gemittelt, gewichtet nach Bevölkerungsgröße, über alle überwachten Kläranlagen.";
      &printGraphicalSynopsis2($path, $text);
    }else{
      print STDERR "Warning: no graphical synopsis stacked timecourse\n";
    }


    if (defined($data{detection}->{Overview}->{Austria})){
      my $path = $data{detection}->{Overview}->{Austria};
      my $text = "Alle im Berichtszeitraum detektierten Varianten, aufgeschlüsselt nach der Anzahl der Kläranlagen mit positivem Signal.";
      &printGraphicalSynopsis2($path, $text);
    }else{
      print STDERR "Warning: no graphical synopsis stacked timecourse\n";
    }

    if (defined($data{sankey}->{Overview}->{Austria})){
      my $path = $data{sankey}->{Overview}->{Austria};
      my $text = "aus dem ganzen Bundesgebiet";
      &printGraphicalSynopsis($path, $text);
    }else{
      print STDERR "Warning: no graphical synopsis on nation level\n";
    }

    if (defined($data{sankey}->{Overview}->{Austria})){
      my $path = $data{sankey}->{Overview}->{allStates};
      my $text = "für jedes Bundesland getrennt";
      &printGraphicalSynopsis($path, $text);
    } else{
      print STDERR "Warning: no graphical synopsis on state level\n";
    }

}
else{
  print STDERR "Warning: no graphical synopsis\n";
}

print $txt3;
}

sub printDMH {
my $txt = '
\newpage
\section{Daten und Methodik}

Die Probenaufbereitung erfolgt bzw. erfolgte durch folgende Kollaborationspartner:
\begin{itemize}
  \item Herbert Oberacher, Institute für Gerichtliche Medizin, Medizinische Universität Innsbruck
  \item Norbert Kreuzinger, Instituts für Wassergüte und Ressourcenmanagement, Technischen Universität Wien
  \item Heribert Insam, Institut für Mikrobiologie, Universität Innsbruck
\end{itemize}

Das Projekt wird bzw. wurde finanziell durch Beiträge von folgenden Partner ermöglicht:
\begin{itemize}
  \item Bundesministerium für Soziales, Gesundheit, Pflege und Konsumentenschutz
  \item Bundesministerium für Bildung, Wissenschaft und Forschung
  \item Agentur für Gesundheit und Ernährungssicherheit
  \item Land Salzburg
  \item Land Vorarlberg
  \item Land Kärnten
\end{itemize}

Ganzgenome Sequenzierung des Sars-CoV-2 Genomes mittels Amplicon-Seq wird am \href{www.biomedical-sequencing.org}{BSF} (Biomedical Sequencing Facility) durchgeführt.

Die Datenanalyse erfolgt am Center for Molecular Medicine (CeMM), in der Arbeitsgruppe Andreas Bergthaler. Die Identifizierung und Quantifizierung der Mutationen im Sars-CoV-2 Genome folgt der etablierten Analyse Pipeline um die Software lofreq (\href{http://doi.org/10.1126/scitranslmed.abe2555 }{Popa A et al. Science Translational Medicine 2020}). \\

Die Methode zur Detektion und Quantifizierung der einzelnen Varianten wurde umfassend validiert und in folgender Publikation beschrieben (\href{https://doi.org/10.1038/s41587-022-01387-y}{Amman, Markt, et al. Nature Biotechnologie 2022}).\\

\subsection{Haftungsausschluss}

Die Analyse Methode ist Gegenstand aktiver Weiterentwicklung. Daher kann, trotz sorgfältiger Analyse, für die Richtigkeit, Vollständigkeit und Aktualität der Angaben keine Garantie übernommen werden.
';
print $txt;
}
