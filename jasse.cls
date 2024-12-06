%\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{jasse}[2014/07/29 ver.0.5]
\LoadClass[a4paper,11pt]{article}
%
\oddsidemargin=27pt
\topmargin=4.6mm
\headheight=10pt
\headsep=10pt
\footskip=30pt
\textheight=600pt
\textwidth=400pt
\columnsep=5mm
%%%%%%%%
\def\@maketitle{%
\newpage\null\vskip 2em%
\begin{center}%
	\let \footnote \thanks%
	{\LARGE\bf \@title \par}%
    	\vskip 1.5em%
    	{\large \lineskip .5em%
    	\begin{tabular}[t]{c}%
		\@author%
	\end{tabular}\par}%
    $^*$\texttt{\@email}\par%
	\vskip 1em
	{\small\@Received; \@Accepted; \@Published}
\end{center}%
  \par%
  \vskip 1.5em%
  \gdef\@email{}%
}%
\long
\def\email#1{\long\gdef\@email{#1}}%
\def\Received#1{\gdef\@Received{Received: #1}}
\def\Accepted#1{\gdef\@Accepted{Accepted: #1}}
\def\Published#1{\gdef\@Published{Published: #1}}
%%%%%%%%%%%%%%%%%%%
% Abstract
%%%%%%%%%%%%%%%%%%%
\def\abstract{%
\noindent{\bfseries{Abstract.}\ }%
\ignorespaces%
}%
%%%%%%%%%%%%%%%%%%%
% Keywords
%%%%%%%%%%%%%%%%%%%
\def\keywords{%
\vskip .5em
\noindent{\bfseries{Keywords:}\ }%
\kern\z@\ignorespaces
}
%%%%%%%%%%%%%%%%%%%
% Acknowledgement
%%%%%%%%%%%%%%%%%%%
\def\acknowledgement{%
\vskip 1.7em
\noindent{\bfseries{\Large Acknowledgement}\ }\\[3mm]
%
}
%%%%%%%%%%%%%%%%%%%
% Appendix
%%%%%%%%%%%%%%%%%%%
\let\oldappendix\appendix
\renewcommand{\appendix}{%
\oldappendix%
\renewcommand{\thesection}{Appendix \Alph{section}.}}%
%%%%%%%%%%%%%%%%%%%
% Include packages
%%%%%%%%%%%%%%%%%%%
%\usepackage{amsmath,amssymb,amsthm}
\usepackage{amsmath,amsthm}
\usepackage[varg]{txfonts}
%\usepackage{amsfonts}
\usepackage{authblk}
\usepackage{lastpage}
\usepackage{fancyhdr}
\usepackage[dvipdfmx]{graphicx}
\renewcommand\Affilfont{\normalsize}% font size of affiliations
\renewcommand\Authands{, }% omit ``, and'' in authors' list 
\renewcommand{\thesection}{\arabic{section}.}
\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}.}
\renewcommand{\thesubsubsection}{\arabic{section}.\arabic{subsection}.\arabic{subsubsection}.}
%
\pagestyle{fancy}
\rhead{\sl Journal of Advanced Simulation in Science and Engineering}
\lhead{}
%
%%%%%%%%%%%%%%%%%%%
% Paper Information
%%%%%%%%%%%%%%%%%%%
\def\Year#1{\gdef\@Year{#1}}
\def\Vol#1{\gdef\@Vol{#1}}
\def\No#1{\gdef\@No{#1}}
\def\Pages#1{\gdef\@Pages{#1}}
\fancypagestyle{titlepage}
   {
   \fancyhead[L]{\footnotesize{\sc J.~Adv.~Simulat.~Sci.~Eng.}\\
   \rm Vol. \@Vol, No. \@No, \@Pages.}%~(\pageref{LastPage}~pages).}
   \fancyhead[R]{\footnotesize\copyright~\@Year\ Japan Society for Simulation Technology}
   \fancyfoot[c]{\thepage}
   \renewcommand\headrulewidth{0pt}
   \renewcommand\footrulewidth{0pt}
   \setlength\headheight{40pt}
   }
%
%%%%%%%%%%%%%%%%%%%
% Theorem environment
%%%%%%%%%%%%%%%%%%%
\newtheorem{theorem}{Theorem}
\newtheorem{algorithm}{Algorithm}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{condition}{Condition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{remark}{Remark}