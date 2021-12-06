
%style file for FastChem manual, based on scrguide from the KOMA script class
\ProvidesClass{fcmanual}[FastChem manual style class]


% Note: You have to pass options to class \textsf{scrguide}, because only
% these options will be processed by |\LoadClassWithOptions| below.
%    \begin{macrocode}
\ExecuteOptions{pagereftolabel,lessindexnotes}
\PassOptionsToClass{%
  bibliography=totoc,
  index=totoc,
  captions=tableabove,
  headings=small,
  headings=optiontoheadandtoc,
  version=last,
}{scrguide}
\ProcessOptions\relax


\LoadClassWithOptions{scrbook}
\RequirePackage[TS1,T1]{fontenc}

\RequirePackage{afterpage}
%\RequirePackage{makeidx}
\RequirePackage{scrdate}
\RequirePackage{scrtime}
%\RequirePackage{picture}
\RequirePackage{graphicx}
%\RequirePackage{booktabs}
%\RequirePackage{longtable}
\RequirePackage{amsmath}
\RequirePackage{listings}
%\RequirePackage{multicol}
%\RequirePackage{marginnote}
\RequirePackage{tabularx}
%\RequirePackage{scrhack}%
%\RequirePackage{mparhack}%


\newif\ifshortversion\shortversionfalse
\InputIfFileExists{scrnonfree.tex}{}{%
  \RequirePackage[svgnames]{xcolor}%
%  \RequirePackage{tocstyle}% hyperref needs xcolor with tocstyle
  \RequirePackage{tikz}%
  \usetikzlibrary{calc}%
  \RequirePackage[open,openlevel=0]{bookmark}% improved pdf bookmarks
 
  \usepackage{lmodern}%
  \@twosidefalse\@openrightfalse\@mparswitchfalse
  \RequirePackage{geometry}%
  \geometry{papersize={180mm,210mm},%
    includehead,includemp,reversemp, marginparwidth=0em, marginparsep=0em,
    vmargin={1.5mm,3mm},hmargin=5.55mm}%
  \RequirePackage{scrlayer-scrpage}%
  \clearpairofpagestyles
  
  \ihead{%
    \smash{%
      \hspace*{.5em}\rule[-\dp\strutbox]{0pt}{\headheight}\headmark
    }%
  }%
  
  \ohead*{\smash{\pagemark\hspace{2mm}}}%
  \colorlet{headbackground}{DarkBlue}%
  \colorlet{headline}{white}%
  \setkomafont{pageheadfoot}{\small\sffamily\bfseries\color{headline}}%
  \setkomafont{pagenumber}{}%
  \DeclareNewLayer[head,background,
    hoffset=0pt,
    width=\paperwidth,
    addheight=\dimexpr \topmargin+1in+\dp\strutbox\relax,
    addvoffset=2\dp\strutbox,
    contents={%
      \color{headbackground}\rule{\layerwidth}{\layerheight}}]
  {headbackground}%
  \AddLayersAtBeginOfPageStyle{scrheadings}{headbackground}%
 
 
  \renewcommand*{\chapterformat}{\thechapter\ }%
  \renewcommand*{\chapterpagestyle}{scrheadings}%
  \renewcommand*{\partpagestyle}{scrheadings}% 
  \setlength{\@tempdimc}{\oddsidemargin}%
  \addtolength{\@tempdimc}{1in}%
  \setheadwidth[-\@tempdimc]{paper}%
  
  \hypersetup{colorlinks,citecolor=Green,linkcolor=Crimson,urlcolor=Blue}%
 
  \let\orig@maketitle\maketitle
  \def\maketitle{%
    \begingroup
      \@twosidetrue
      \let\cleardoublepage\clearpage
      \addtolength{\oddsidemargin}{\evensidemargin}%
      \setlength{\evensidemargin}{.5\oddsidemargin}%
      \setlength{\oddsidemargin}{\evensidemargin}%
      \orig@maketitle\relax
    \endgroup
  }%
  \renewcommand*{\chapterheadstartvskip}{}%
  \raggedbottom
  \sloppy
  \lstset{breaklines,prebreak=\mbox{$\hookleftarrow$}}%
}

\RequirePackage{textcomp}


\newcommand*{\indexrm}[1]{\textrm{\hyperpage{#1}}}
\newcommand*{\indexother}{}\let\indexother\indexrm
\newcommand*{\indexit}[1]{\textit{\hyperpage{#1}}}
\newcommand*{\indexbf}[1]{\textbf{\hyperpage{#1}}}
\newcommand*{\indexmain}{}\let\indexmain\indexbf
\newcommand*{\indexsl}[1]{\textsl{\hyperpage{#1}}}
\newcommand*{\indexsf}[1]{\textsf{\hyperpage{#1}}}
\newcommand*{\indexsc}[1]{\textsc{\hyperpage{#1}}}


\let\jnl@style=\rmfamily 
\def\ref@jnl#1{{\jnl@style#1}}% 
\newcommand\aj{\ref@jnl{AJ}}%        % Astronomical Journal 
\newcommand\araa{\ref@jnl{ARA\&A}}%  % Annual Review of Astron and Astrophys 
\newcommand\apj{\ref@jnl{ApJ}}%    % Astrophysical Journal ++
\newcommand\apjl{\ref@jnl{ApJL}}     % Astrophysical Journal, Letters 
\newcommand\apjs{\ref@jnl{ApJS}}%    % Astrophysical Journal, Supplement 
\newcommand\ao{\ref@jnl{ApOpt}}%   % Applied Optics ++
\newcommand\apss{\ref@jnl{Ap\&SS}}%  % Astrophysics and Space Science 
\newcommand\aap{\ref@jnl{A\&A}}%     % Astronomy and Astrophysics 
\newcommand\aapr{\ref@jnl{A\&A~Rv}}%  % Astronomy and Astrophysics Reviews 
\newcommand\aaps{\ref@jnl{A\&AS}}%    % Astronomy and Astrophysics, Supplement 
\newcommand\azh{\ref@jnl{AZh}}%       % Astronomicheskii Zhurnal 
\newcommand\baas{\ref@jnl{BAAS}}%     % Bulletin of the AAS 
\newcommand\icarus{\ref@jnl{Icarus}}% % Icarus
\newcommand\jrasc{\ref@jnl{JRASC}}%   % Journal of the RAS of Canada 
\newcommand\memras{\ref@jnl{MmRAS}}%  % Memoirs of the RAS 
\newcommand\mnras{\ref@jnl{MNRAS}}%   % Monthly Notices of the RAS 
\newcommand\pra{\ref@jnl{PhRvA}}% % Physical Review A: General Physics ++
\newcommand\prb{\ref@jnl{PhRvB}}% % Physical Review B: Solid State ++
\newcommand\prc{\ref@jnl{PhRvC}}% % Physical Review C ++
\newcommand\prd{\ref@jnl{PhRvD}}% % Physical Review D ++
\newcommand\pre{\ref@jnl{PhRvE}}% % Physical Review E ++
\newcommand\prl{\ref@jnl{PhRvL}}% % Physical Review Letters 
\newcommand\pasp{\ref@jnl{PASP}}%     % Publications of the ASP 
\newcommand\pasj{\ref@jnl{PASJ}}%     % Publications of the ASJ 
\newcommand\qjras{\ref@jnl{QJRAS}}%   % Quarterly Journal of the RAS 
\newcommand\skytel{\ref@jnl{S\&T}}%   % Sky and Telescope 
\newcommand\solphys{\ref@jnl{SoPh}}% % Solar Physics 
\newcommand\sovast{\ref@jnl{Soviet~Ast.}}% % Soviet Astronomy 
\newcommand\ssr{\ref@jnl{SSRv}}% % Space Science Reviews 
\newcommand\zap{\ref@jnl{ZA}}%       % Zeitschrift fuer Astrophysik 
\newcommand\nat{\ref@jnl{Nature}}%  % Nature 
\newcommand\iaucirc{\ref@jnl{IAUC}}% % IAU Cirulars 
\newcommand\aplett{\ref@jnl{Astrophys.~Lett.}}%  % Astrophysics Letters 
\newcommand\apspr{\ref@jnl{Astrophys.~Space~Phys.~Res.}}% % Astrophysics Space Physics Research 
\newcommand\bain{\ref@jnl{BAN}}% % Bulletin Astronomical Institute of the Netherlands 
\newcommand\fcp{\ref@jnl{FCPh}}%   % Fundamental Cosmic Physics 
\newcommand\gca{\ref@jnl{GeoCoA}}% % Geochimica Cosmochimica Acta 
\newcommand\grl{\ref@jnl{Geophys.~Res.~Lett.}}%  % Geophysics Research Letters 
\newcommand\jcp{\ref@jnl{JChPh}}%     % Journal of Chemical Physics 
\newcommand\jgr{\ref@jnl{J.~Geophys.~Res.}}%     % Journal of Geophysics Research 
\newcommand\jqsrt{\ref@jnl{JQSRT}}%   % Journal of Quantitiative Spectroscopy and Radiative Trasfer 
\newcommand\memsai{\ref@jnl{MmSAI}}% % Mem. Societa Astronomica Italiana 
\newcommand\nphysa{\ref@jnl{NuPhA}}%     % Nuclear Physics A 
\newcommand\physrep{\ref@jnl{PhR}}%       % Physics Reports 
\newcommand\physscr{\ref@jnl{PhyS}}%        % Physica Scripta 
\newcommand\planss{\ref@jnl{Planet.~Space~Sci.}}%  % Planetary Space Science 
\newcommand\procspie{\ref@jnl{Proc.~SPIE}}%      % Proceedings of the SPIE 

\newcommand\actaa{\ref@jnl{AcA}}%  % Acta Astronomica
\newcommand\caa{\ref@jnl{ChA\&A}}%  % Chinese Astronomy and Astrophysics
\newcommand\cjaa{\ref@jnl{ChJA\&A}}%  % Chinese Journal of Astronomy and Astrophysics
\newcommand\jcap{\ref@jnl{JCAP}}%  % Journal of Cosmology and Astroparticle Physics
\newcommand\na{\ref@jnl{NewA}}%  % New Astronomy
\newcommand\nar{\ref@jnl{NewAR}}%  % New Astronomy Review
\newcommand\pasa{\ref@jnl{PASA}}%  % Publications of the Astron. Soc. of Australia
\newcommand\rmxaa{\ref@jnl{RMxAA}}%  % Revista Mexicana de Astronomia y Astrofisica

%% added feb 9, 2016
\newcommand\maps{\ref@jnl{M\&PS}}% Meteoritics and Planetary Science
\newcommand\aas{\ref@jnl{AAS Meeting Abstracts}}% American Astronomical Society Meeting Abstracts
\newcommand\dps{\ref@jnl{AAS/DPS Meeting Abstracts}}% American Astronomical Society/Division for Planetary Sciences Meeting Abstracts



\let\astap=\aap 
\let\apjlett=\apjl 
\let\apjsupp=\apjs 
\let\applopt=\ao 

\endinput
