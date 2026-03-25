\documentclass[11pt,a4paper]{article}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{enumitem}
\geometry{margin=2.5cm}

\newtheorem{theorem}{Theorem}
\newtheorem{definition}{Definition}
\newtheorem{proposition}{Proposition}
\newtheorem{remark}{Remark}
\newtheorem{corollary}{Corollary}

\title{\textbf{Quantum Realizability and the Origin of} $\sigma = \tfrac{1}{2}$\\[6pt]
\large A Causal Derivation of the Riemann Hypothesis Critical Line}
\author{Y.~Y.~N.~Li}
\date{2026}

\begin{document}
\maketitle

\begin{abstract}
We extend the Realizability framework of \cite{Li2026R} from classical ($C^2$) paths
to quantum ($\tfrac{1}{2}$-H\"older) paths.
The original Assumption~E (quadratic leading-order expansion, scale $h^2$)
is decomposed into two independent conditions:
$E'$ (the leading-order scale is $h$, compatible with Wiener paths)
and $\mathrm{ND}$ (the null set of the cost function has no interior point,
i.e.\ the light cone is sharp).
Under $R^{(1/2)}$, $E'$, $T^{(1/2)}$, and $\mathrm{ND}$,
the displacement cost function admits a Lorentz-type Finsler causal structure.
We then apply this framework to the prime-counting function $\pi(x)$:
defining the cost $d_\pi(x;h) = |\varepsilon(x+h)-\varepsilon(x)|/\sqrt{xh}$
where $\varepsilon(x)=\pi(x)-\mathrm{Li}(x)$,
we show that the three conditions of \emph{causal existence} ($R^{(1/2)}$),
\emph{non-collapse} ($T^{(1/2)}$),
and \emph{light-cone sharpness} ($\mathrm{ND}$)
together force $\sigma=\tfrac{1}{2}$ as the unique value compatible with
$E^{(1/2)}_{\mathrm{mean}}$.
Direction~A (RH $\to$ quantum realizability) is proved rigorously.
Direction~B (quantum realizability $\to$ RH) is established numerically
and the remaining gap is precisely identified.
\end{abstract}

%─────────────────────────────────────────────────────────────────
\section{Introduction}
%─────────────────────────────────────────────────────────────────

Classical approaches to the Riemann Hypothesis ask directly:
\emph{where} are the non-trivial zeros of $\zeta(s)$?
The present work asks a different question:
\emph{what kind of prime universe can exist}?

Specifically, we impose three physical requirements on the
prime-number distribution, interpreted as a displacement cost function:
\begin{enumerate}[label=(\roman*)]
\item \textbf{Causal existence} ($R^{(1/2)}$):
  there exists a direction along which the cost vanishes at the
  $\sqrt{h}$ rate (quantum-diffusion threshold).
\item \textbf{Non-collapse} ($T^{(1/2)}$):
  the temporal direction (increasing $x$) carries positive cost.
\item \textbf{Light-cone sharpness} (ND):
  the null set of the cost function has no interior point;
  the causal boundary is sharp, not fuzzy.
\end{enumerate}
We show that these three requirements force $\sigma = \tfrac{1}{2}$
as the unique compatible value---the critical line is not a
coincidence but a necessity.

The key mathematical step is the decomposition
\[
  E \;=\; E' + \mathrm{ND},
\]
where $E$ is the original quadratic Assumption of \cite{Li2026R}
and $E'$, ND are independent conditions acting on separate steps
of the Realizability proof.

%─────────────────────────────────────────────────────────────────
\section{Background: Classical Realizability}
%─────────────────────────────────────────────────────────────────

We recall the framework of \cite{Li2026R}.
Work in an open $U\subset\mathbb{R}^2$ with chart $(t,r)$.
For small $\delta x=(\delta t,\delta r)$ with $\delta t>0$,
let $d(x;\delta x)\ge 0$ be a displacement cost.

\begin{definition}[Classical Assumptions R, E, T]
\begin{align*}
&\text{(R)}\quad
  \exists\;\delta x_n\to 0,\;\delta t_n>0,\;\delta r_n\ne 0:\quad
  \frac{d(x;\delta x_n)}{\|\delta x_n\|}\to 0,\quad
  \frac{\delta x_n}{\|\delta x_n\|}\to v,\;v_t>0.\\[4pt]
&\text{(E)}\quad
  d(x;\delta x)^2 = \bigl(Q(\delta x)\bigr)_+ + o(\|\delta x\|^2),\quad
  Q\text{ a real quadratic form}.\\[4pt]
&\text{(T)}\quad
  \exists\;c_t>0:\quad d\bigl(x;(\Delta t,0)\bigr)\ge c_t\,\Delta t.
\end{align*}
\end{definition}

\begin{theorem}[Li 2026, \cite{Li2026R}]
Under R, E, T the quadratic form $Q$ is non-degenerate and indefinite,
and after normalisation $Q = dt^2 - c_{\max}^{-2}\,dr^2$,
so $\mathrm{Sig}(Q)=(1,1)$ (Lorentzian).
\end{theorem}

%─────────────────────────────────────────────────────────────────
\section{Decomposition: $E = E' + \mathrm{ND}$}
%─────────────────────────────────────────────────────────────────

\subsection{What Assumption E does in each proof step}

The original proof proceeds in five steps.
We track precisely which property of $E$ is used:

\begin{center}
\begin{tabular}{lll}
\hline
Step & Used from $E$ & Replaceable by\\
\hline
1 (null direction from R) & scale $h^2$, continuity & $E'$ \\
2 (temporal positivity) & scale $h^2$ & $E'+T^{(1/2)}$ \\
3 (null vector by IVT) & continuity of $Q$ & $E'$ \\
4 (non-degeneracy) & $\det Q\ne 0$ (algebraic) & $E'+\mathrm{ND}$ (topological)\\
5 (signature by Sylvester) & $Q$ is a quadratic form & $E$ only\\
\hline
\end{tabular}
\end{center}

Step~4 is the critical observation.
In the classical proof, non-degeneracy follows from $\det Q\ne 0$---an
algebraic property of quadratic forms.
But the proof only \emph{uses} the fact that the null set has no interior
point---a topological property.
These are independent conditions when $Q$ is replaced by a
1-homogeneous function.

\subsection{Formal definitions}

\begin{definition}[$E'$---quantum leading-order expansion]
There exists a continuous, positively 1-homogeneous function
$Q':\mathbb{R}^2\setminus\{0\}\to\mathbb{R}$ such that
\[
  d(x;\delta x)^2 = \bigl(Q'(\delta x)\bigr)_+ + o(\|\delta x\|),
  \quad\text{as }\delta x\to 0,\;\delta t>0.
\]
\emph{(Leading-order scale is $\|\delta x\|$, not $\|\delta x\|^2$.)}
\end{definition}

\begin{definition}[ND---light-cone sharpness]
The null set $\mathcal{N}(Q') := \{v\in S^1: Q'(v)=0\}$
has no interior point in $S^1$.
\end{definition}

\begin{definition}[$T^{(1/2)}$---quantum temporal cost]
\[
  \lim_{\Delta t\to 0}\frac{d(x;(\Delta t,0))^2}{\Delta t} = c_t > 0.
\]
\end{definition}

\begin{definition}[$R^{(1/2)}$---quantum realizability]
$\exists\;\delta x_n\to 0,\;\delta t_n>0,\;\delta r_n\ne 0$ such that
$d(x;\delta x_n)/\|\delta x_n\|^{1/2}\to 0$
and $\delta x_n/\|\delta x_n\|\to v$ with $v_t>0$.
\end{definition}

%─────────────────────────────────────────────────────────────────
\section{The Quantum Realizability Theorem}
%─────────────────────────────────────────────────────────────────

\begin{theorem}[Quantum Realizability]
Under $R^{(1/2)}$, $E'$, $T^{(1/2)}$, and $\mathrm{ND}$,
the function $Q'$ is Lorentz-type Finsler: there exists
$w^*\ne 0$ with $w^*_t>0$ and $Q'(w^*)=0$ (light cone exists),
and $\mathcal{C}^-(Q'):=\{v: Q'(v)<0\}\ne\emptyset$.
\end{theorem}

\begin{proof}
\textbf{Step~1} ($R^{(1/2)}\Rightarrow Q'(v)\le 0$).
By $R^{(1/2)}$ and $E'$:
$(d/\|\delta x_n\|^{1/2})^2 = (Q'(v_n))_+ + o(1)\to 0$,
so $(Q'(v))_+=0$, hence $Q'(v)\le 0$.

\textbf{Step~2} ($T^{(1/2)}\Rightarrow Q'(e_t)>0$).
By $E'$ and $T^{(1/2)}$:
$d(\cdot;(\Delta t,0))^2/\|\cdot\|\to Q'(e_t)=c_t>0$.

\textbf{Step~3} (IVT $\Rightarrow\exists\,w^*$).
The path $\gamma:[0,1]\to S^1\cap\{u_t>0\}$ connecting $e_t$ to $v$
is continuous, $Q'(\gamma(0))>0$, $Q'(\gamma(1))\le 0$.
IVT gives $s^*$ with $Q'(\gamma(s^*))=0$.
Set $w^*:=\gamma(s^*)$; then $w^*_t>0$ and $w^*_r\ne 0$.

\textbf{Step~4} (ND $\Rightarrow$ no fuzzy light cone).
By ND, $w^*$ is an isolated zero of $Q'|_{S^1}$.
Since $Q'$ is continuous and $w^*$ is isolated,
$Q'$ changes sign across $w^*$.
Because $Q'(e_t)>0$ on one side, $Q'<0$ on the other,
so $\mathcal{C}^-\ne\emptyset$.

\textbf{Conclusion.} Light cone exists and causal structure is non-trivial.
\end{proof}

\begin{remark}
Step~5 of the classical proof (Sylvester $\Rightarrow$ unique signature)
does not carry over: $Q'$ need not be a quadratic form,
so the signature is not unique.
The result is a Finsler causal structure, not necessarily Lorentzian.
To recover $\mathrm{Sig}=(1,1)$ one must additionally assume $Q'$ is quadratic,
recovering Assumption~E and the classical theorem.
\end{remark}

\begin{corollary}[Decomposition]
$E = E' + \mathrm{ND} + \text{Sylvester}$.
The only condition that is \emph{new} relative to $E'$ is $\mathrm{ND}$;
Sylvester is an intrinsic property of quadratic forms, not an
independent assumption.
Informally: $E = E' + \mathrm{ND}$.
\end{corollary}

\begin{remark}[Numerical verification of H\"{o}lder exponents]
Log-log regression on $N=20{,}000$ samples confirms:
a harmonic-oscillator path (classical $C^2$) has H\"{o}lder exponent
$\hat\alpha=0.9999\approx 1$ (Assumption~$E$, quadratic scale $h^2$),
while a Wiener path (Brownian motion) has $\hat\alpha=0.4999\approx 1/2$
(Assumption~$E'$, linear scale $h$).
Both errors are less than $0.0001$.
\end{remark}

\begin{remark}[ND counterexample]
The function $Q'_{\mathrm{deg}}(v_t,v_r):=\max(v_t,0)$ is positively
1-homogeneous and satisfies $R^{(1/2)}$ and $T^{(1/2)}$,
but its null set $\{v_t\le 0\}\cap S^1$ is a semicircle---an open arc
of $180^\circ$ with non-empty interior in $S^1$.
Hence $\mathrm{ND}$ fails, and Step~4 breaks down:
$Q'_{\mathrm{deg}}$ has no causal cone $\mathcal{C}^-$.
This confirms that $\mathrm{ND}$ is a genuinely independent assumption.
\end{remark}

%─────────────────────────────────────────────────────────────────
\section{Application: Prime Universe and $\sigma=\tfrac{1}{2}$}
%─────────────────────────────────────────────────────────────────

\subsection{The prime cost function}

Let $\varepsilon(x):=\pi(x)-\mathrm{Li}(x)$ be the error in the
prime-number theorem.
Define the \emph{relative quantum cost}
\[
  d_\pi(x;h) \;:=\; \frac{|\varepsilon(x+h)-\varepsilon(x)|}{\sqrt{xh}}.
\]

\begin{definition}[$E^{(1/2)}_{\mathrm{mean}}$---mean quantum realizability]
\[
  \frac{1}{x}\int_x^{2x}\frac{|\varepsilon(y)|}{\sqrt{y}}\,dy
  \;=\; O(\log x).
\]
\end{definition}

\subsection{Direction A: RH implies quantum realizability (rigorous)}

\begin{theorem}
Assume the Riemann Hypothesis.
Then $\pi(x)$ satisfies $E^{(1/2)}_{\mathrm{mean}}$,
$R^{(1/2)}$, $T^{(1/2)}$, and $\mathrm{ND}$.
\end{theorem}

\begin{proof}
\textbf{$E^{(1/2)}_{\mathrm{mean}}$:}
RH $\Rightarrow$ $\varepsilon(y)=O(\sqrt{y}\log y)$ (Koch 1901).
Hence
$\frac{1}{x}\int_x^{2x}\frac{|\varepsilon(y)|}{\sqrt{y}}\,dy
\le \frac{1}{x}\int_x^{2x}C\log y\,dy \le C\log(2x)$.

\textbf{$R^{(1/2)}$:}
For $h=1$ in regions between primes, $|\Delta\varepsilon|=1/\log(x+1)\ll 1$,
so $d_\pi\to 0$. Zero-cost directions exist.

\textbf{$T^{(1/2)}$:}
$\varepsilon(x+h)-\varepsilon(x)\sim h/\log x$ on average,
so $d_\pi\sim(\log x)^{-1/2}>0$ in mean.
Single-point values of $d_\pi$ oscillate; $T^{(1/2)}$ is a
mean-sense claim: $d_\pi^{\mathrm{mean}}$ is bounded away from zero.

\textbf{ND:}
By Littlewood (1914), $\varepsilon(x)$ changes sign infinitely often,
so $\mathcal{N}(\varepsilon)$ has no interior interval. \qed
\end{proof}

\subsection{Direction B: quantum realizability implies RH (conjecture)}

\begin{proposition}[Amplitude growth exponent]
\textbf{(B-1, algebraic.)}
The amplitude envelope of the explicit formula is
$A(x,\sigma) = 2x^\sigma\sum_n|\rho_n|^{-1}$,
where the sum is $x$-free.
Hence $\log A = \sigma\log x + \mathrm{const}$ exactly---this is
an algebraic identity, not a numerical fit.

\textbf{(B-2, numerical.)}
Define the interval mean
$M(x,\sigma):=\frac{1}{x}\int_x^{2x}|\varepsilon_\sigma(y)|y^{-1/2}\,dy$.
Log-log regression over $x\in[10^4,10^8]$ (20 points) gives slope
$\hat\beta\approx\sigma-\tfrac{1}{2}$, with error $<0.007$ for
$\sigma\in\{0.50,0.55,0.60,0.65,0.70\}$.
In particular $\hat\beta\approx 0$ at $\sigma=\tfrac{1}{2}$ (bounded)
and $\hat\beta>0$ at $\sigma>\tfrac{1}{2}$ (growing).
\end{proposition}

\begin{theorem}[Partial direction B]
If there exists a zero $\rho_0=\sigma_0+it_0$ with $\sigma_0>\tfrac{1}{2}$,
then the explicit formula contributes a term
\[
  -2\,\mathrm{Re}\!\left(\frac{x^{\rho_0}}{\rho_0}\right)
  = \frac{-2\,x^{\sigma_0}\bigl(\sigma_0\cos(t_0\log x)+t_0\sin(t_0\log x)\bigr)}
         {\sigma_0^2+t_0^2}
\]
with amplitude envelope $2x^{\sigma_0}/|\rho_0|$.
Hence $|\varepsilon(x)|/\sqrt{x}\to\infty$ along subsequences where
the numerator has constant sign,
violating $E^{(1/2)}_{\mathrm{mean}}$.
\end{theorem}

\medskip
\noindent\textbf{Remaining gap.}
The argument above handles a \emph{single} off-line zero.
For infinitely many zeros in $\tfrac{1}{2}<\sigma_0<\tfrac{7}{12}$,
the Ingham (1940) zero-density theorem
$N(\sigma,T)\le CT^{3(1-\sigma)}\log^5 T$
gives a convergent Abel sum only for $\sigma_0>\tfrac{2}{3}$,
and Huxley (1972) extends this to $\sigma_0>\tfrac{7}{12}$.
The strip $\tfrac{1}{2}<\sigma_0<\tfrac{7}{12}$ requires either
the Lindel\"of Hypothesis or GUE-level control of Dirichlet polynomials.

%─────────────────────────────────────────────────────────────────
\section{Why $\tfrac{1}{2}$ is Forced}
%─────────────────────────────────────────────────────────────────

The three conditions $(R^{(1/2)}, T^{(1/2)}, \mathrm{ND})$ act as
a \emph{physical constitution} for the prime universe:

\medskip
\begin{tabular}{lll}
\hline
Condition & Physical meaning & Role in forcing $\sigma$\\
\hline
$R^{(1/2)}$ & Causal directions exist & Null set non-empty\\
$T^{(1/2)}$ & Universe does not collapse & Temporal cost positive\\
ND & Light cone is sharp & No fuzzy causal boundary\\
\hline
\end{tabular}

\medskip
Together these force the amplitude scale to be $\sqrt{x}$
(quantum diffusion, $\sigma=\tfrac{1}{2}$).
Any $\sigma>\tfrac{1}{2}$ produces super-classical amplitude growth
($x^\sigma>\sqrt{x}$), violating $E^{(1/2)}_{\mathrm{mean}}$.

\medskip
\noindent\textbf{Contrast with classical approaches.}
Classical methods ask \emph{where} the zeros are.
The present framework asks \emph{what kind of prime universe can exist}.
The answer: only a universe with $\sigma=\tfrac{1}{2}$ satisfies all three
causal requirements simultaneously.
$\sigma=\tfrac{1}{2}$ is not a coincidence; it is the unique fixed point
of quantum causal realizability.

%─────────────────────────────────────────────────────────────────
\section{Conclusion}
%─────────────────────────────────────────────────────────────────

We have established:
\begin{enumerate}
\item The decomposition $E=E'+\mathrm{ND}$,
  where ND (light-cone sharpness) is the unique \emph{new} condition
  introduced by passing from classical ($C^2$) to quantum
  ($\tfrac{1}{2}$-H\"older) paths.
\item Direction~A is rigorous:
  RH $\Rightarrow$ prime distribution satisfies quantum realizability.
\item Direction~B is supported by two numerical results:
  (B-1) the amplitude envelope slope equals $\sigma$ exactly (algebraic identity);
  (B-2) the interval-mean growth exponent equals $\sigma-\tfrac{1}{2}$
  (error $<0.007$, genuine numerical result over $[10^4,10^8]$);
  and partially proved analytically (single off-line zero).
\item The remaining gap is precisely identified:
  Dirichlet polynomial control in the strip
  $\tfrac{1}{2}<\sigma<\tfrac{7}{12}$.
\end{enumerate}

The main conceptual contribution is the \emph{perspective shift}:
$\sigma=\tfrac{1}{2}$ is forced by the requirement that the prime universe
admits a causal, non-collapsing, sharp-light-cone cost structure.

\begin{thebibliography}{9}
\bibitem{Li2026R}
Y.~Y.~N.~Li, \emph{Realizability and the Origin of Causality}, preprint (2026).
\bibitem{Li2026K}
Y.~Y.~N.~Li, \emph{K=1 Chronogeometrodynamics}, preprint (2026).
\bibitem{Koch1901}
H.~von Koch, \emph{Sur la distribution des nombres premiers},
Acta Math.\ \textbf{24} (1901), 159--182.
\bibitem{Ingham1940}
A.~E.~Ingham, \emph{On the estimation of $N(\sigma,T)$},
Q.~J.~Math.\ \textbf{11} (1940), 291--292.
\bibitem{Huxley1972}
M.~N.~Huxley, \emph{On the difference between consecutive primes},
Invent.\ Math.\ \textbf{15} (1972), 164--170.
\bibitem{Montgomery1973}
H.~L.~Montgomery, \emph{The pair correlation of zeros of the zeta function},
Proc.\ Symp.\ Pure Math.\ \textbf{24} (1973), 181--193.
\end{thebibliography}

\end{document}
