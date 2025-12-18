# The theory behind onlineFDR

## FDR Control

Consider a sequence of hypotheses $H_{1},H_{2},H_{3},\ldots$ that arrive
sequentially in a stream, with corresponding $p$-values
$\left( p_{1},p_{2},p_{3},\ldots \right)$. A testing procedure provides
a sequence of adjusted significance thresholds $\alpha_{i}$, with
corresponding decision rule: $$R_{i} = \begin{Bmatrix}
1 & {{\text{if}\mspace{6mu}}p_{i} \leq \alpha_{i}} & \left( {\text{reject}\mspace{6mu}}H_{i} \right) \\
0 & {\text{otherwise}\mspace{6mu}} & \left( {\text{accept}\mspace{6mu}}H_{i} \right)
\end{Bmatrix}$$

In *online* testing, the significance thresholds can only be functions
of the prior decisions,
i.e. $\alpha_{i} = \alpha_{i}\left( R_{1},R_{2},\ldots,R_{i - 1} \right)$.

Javanmard and Montanari (2015, 2018) proposed two procedures for online
control. The first is LOND, which stands for (significance) Levels based
On Number of Discoveries. The second is LORD, which stands for
(significance) Levels based On Recent Discovery. LORD was subsequently
extended by Ramdas *et al.* (2017). Ramdas *et al.* (2018) also proposed
the SAFFRON procedure, which provides an adaptive method of online FDR
control, which includes a variant of Alpha-investing. Finally, Tian &
Ramdas (2019) proposed the ADDIS procedure as an improvement of SAFFRON
in the presence of conservative nulls.

### LOND

The LOND procedure controls the FDR for independent or positively
dependent (PRDS) $p$-values. Given an overall significance level
$\alpha$, we choose a sequence of non-negative numbers
$\beta = \left( \beta_{i} \right)_{i \in {\mathbb{N}}}$ such that they
sum to $\alpha$. The values of the adjusted significance thresholds
$\alpha_{i}$ are chosen as follows:
$$\alpha_{i} = \beta_{i}\left( D(i - 1) + 1 \right)$$ where
$D(n) = \sum_{i = 1}^{n}R_{i}$ denotes the number of discoveries
(i.e. rejections) in the first $n$ hypotheses tested.

LOND can be adjusted to also control FDR under arbitrarily dependent
$p$-values. To do so, it is modified with
${\widetilde{\beta}}_{i} = \beta_{i}/H(i)$ in place of $\beta_{i}$,
where $H(i) = \sum_{j = 1}^{i}\frac{1}{j}$ is the $i$-th harmonic
number. Note that this leads to a substantial loss in power compared to
the unadjusted LOND procedure. The correction factor is similar to the
classical one used by Benjamini and Yekutieli (2001), except that in
this case the $i$-th hypothesis among $N$ is penalised by a factor of
$H(i)$ to give consistent results across time (as compared to a factor
$H(N)$ for the Benjamini and Yekutieli method).

The default sequence of $\beta$ is given by
$$\beta_{j} = C\alpha\frac{\log\left( \max(j,2) \right)}{je^{\sqrt{\log j}}}$$
where $C \approx 0.07720838$, as proposed by Javanmard and Montanari
(2018) equation 31.

### LORD

The LORD procedure controls the FDR for independent $p$-values. We first
fix a sequence of non-negative numbers
$\gamma = \left( \gamma_{i} \right)_{i \in {\mathbb{N}}}$ such that
$\gamma_{i} \geq \gamma_{j}$ for $i \leq j$ and
$\sum_{i = 1}^{\infty}\gamma_{i} = 1$. At each time $i$, let $\tau_{i}$
be the last time a discovery was made before $i$:
$$\tau_{i} = \max\left\{ l \in \{ 1,\ldots,i - 1\}:R_{l} = 1 \right\}$$

LORD depends on constants $w_{0}$ and $b_{0}$, where $w_{0} \geq 0$
represents the initial ‘wealth’ of the procedure and $b_{0} > 0$
represents the ‘payout’ for rejecting a hypothesis. We require
$w_{0} + b_{0} \leq \alpha$ for FDR control to hold.

Javanmard and Montanari (2018) presented three different versions of
LORD, which have different definitions of the adjusted significance
thresholds $\alpha_{i}$. Versions 1 and 2 have since been superseded by
the LORD++ procedure of Ramdas *et al.* (2017), so we do not describe
them here.

- **LORD++**: The significance thresholds for LORD++ are chosen as
  follows:
  $$\alpha_{i} = \gamma_{i}w_{0} + \left( \alpha - w_{0} \right)\gamma_{i - \tau_{1}} + \alpha\sum\limits_{j:\tau_{j} < i,\tau_{j} \neq \tau_{1}}\gamma_{i - \tau_{j}}$$

- **LORD 3**: The significance thresholds depend on the time of the last
  discovery time and the wealth accumulated at that time, with
  $$\alpha_{i} = \gamma_{i - \tau_{i}}W\left( \tau_{i} \right)$$ where
  $\tau_{1} = 0$. Here $\{ W(j)\}_{j \geq 0}$ represents the ‘wealth’
  available at time $j$, and is defined recursively:

$$\begin{aligned}
{W(0)} & {= w_{0}} \\
{W(j)} & {= W(j - 1) - \alpha_{j - 1} + b_{0}R_{j}}
\end{aligned}$$

- **D-LORD**: This is equivalent to the LORD++ procedure with
  discarding. Given a discarding threshold $\tau \in (0,1)$ and initial
  wealth $w_{0} \leq \tau\alpha$ the significance thresholds are chosen
  as follows: $$\alpha_{t} = \min\{\tau,{\widetilde{\alpha}}_{t}\}$$
  where
  $${\widetilde{\alpha}}_{t} = w_{0}\gamma_{S^{t}} + \left( \tau\alpha - w_{0} \right)\gamma_{S^{t} - \kappa_{1}^{*}} + \tau\alpha\sum\limits_{j \geq 2}\gamma_{S^{t} - \kappa_{j}^{*}}$$
  and
  $$\kappa_{j} = \min\{ i \in \lbrack t - 1\rbrack:\sum\limits_{k \leq i}1\{ p_{k} \leq \alpha_{k}\} \geq j\},\;\kappa_{j}^{*} = \sum\limits_{i \leq \kappa_{j}}1\{ p_{i} \leq \tau\},\; S^{t} = \sum\limits_{i < t}1\{ p_{i} \leq \tau\}$$

LORD++ is an instance of a monotone rule, and provably controls the FDR
for independent p-values provided $w_{0} \leq \alpha$. LORD 3 is a
non-monotone rule, and FDR control is only demonstrated empirically. In
some scenarios with large $N$, LORD 3 will have a slightly higher power
than LORD++ (see Robertson *et al.*, 2018), but since it is a
non-monotone rule we would recommend using LORD++ (which is the
default), especially since it also has a provable guarantee of FDR
control.

In all versions, the default sequence of $\gamma$ is given by
$$\gamma_{j} = C\frac{\log\left( \max(j,2) \right)}{je^{\sqrt{\log j}}}$$
where $C \approx 0.07720838$, as proposed by Javanmard and Montanari
(2018) equation 31.

Javanmard and Montanari (2018) also proposed an adjusted version of LORD
that is valid for arbitrarily *dependent* p-values. Similarly to LORD 3,
the adjusted significance thresholds are set equal to
$$\alpha_{i} = \xi_{i}W\left( \tau_{i} \right)$$ where (assuming
$w_{0} \leq b_{0}$),
$\sum_{j = 1}^{\infty}\xi_{i}\left( 1 + \log(j) \right) \leq \alpha/b_{0}$

The default sequence of $\xi$ is given by
$$\xi_{j} = \frac{C\alpha}{b_{0}j\log\left( \max(j,2) \right)^{3}}$$
where $C \approx 0.139307$.

Note that allowing for dependent p-values can lead to a substantial loss
in power compared with the LORD procedures described above.

### SAFFRON

The SAFFRON procedure controls the FDR for independent p-values, and was
proposed by Ramdas *et al.* (2018). The algorithm is based on an
estimate of the proportion of true null hypotheses. More precisely,
SAFFRON sets the adjusted test levels based on an estimate of the amount
of alpha-wealth that is allocated to testing the true null hypotheses.

SAFFRON depends on constants $w_{0}$ and $\lambda$, where $w_{0}$
satisfies $0 \leq w_{0} \leq \alpha$ and represents the initial ‘wealth’
of the procedure, and $\lambda \in (0,1)$ represents the threshold for a
‘candidate’ hypothesis. A ‘candidate’ refers to p-values smaller than
$\lambda$, since SAFFRON will never reject a p-value larger than
$\lambda$. These candidates can be thought of as the hypotheses that are
a-priori more likely to be non-null.

The SAFFRON procedure runs as follows:

1.  At each time $t$, define the number of candidates after the $j$-th
    rejection as
    $$C_{j +} = C_{j +}(t) = \sum\limits_{i = \tau_{j} + 1}^{t - 1}C_{i}$$
    where $C_{t} = 1\{ p_{t} \leq \lambda\}$ is the indicator for
    candidacy.

2.  SAFFRON starts with
    $\alpha_{1} = \min\{(1 - \lambda)\gamma_{1}w_{0},\lambda\}$.
    Subsequent test levels are chosen as
    $\alpha_{t} = \min\{\lambda,{\widetilde{\alpha}}_{t}\}$, where
    $${\widetilde{\alpha}}_{t} = (1 - \lambda)\left\lbrack w_{0}\gamma_{t - C_{0 +}} + \left( \alpha - w_{0} \right)\gamma_{t - \tau_{1} - C_{1 +}} + \alpha\sum\limits_{j \geq 2}\gamma_{t - \tau_{j} - C_{j +}} \right\rbrack$$

The default sequence of $\gamma$ for SAFFRON is given by
$\gamma_{j} \propto j^{- 1.6}$.

### Alpha-investing

Ramdas et al. (2018) proposed a variant of the Alpha-investing algorithm
of Foster and Stine (2008) that guarantees FDR control for independent
p-values. This procedure uses SAFFRON’s update rule with the constant
$\lambda$ replaced by a sequence $\lambda_{i} = \alpha_{i}$. This is
also equivalent to using the ADDIS algorithm (see below) with $\tau = 1$
and $\lambda_{i} = \alpha_{i}$.

### ADDIS

The ADDIS procedure controls the FDR for independent p-values, and was
proposed by Tian & Ramdas (2019). The algorithm compensates for the
power loss of SAFFRON with conservative nulls, by including both
adaptivity in the fraction of null hypotheses (like SAFFRON) and the
conservativeness of nulls (unlike SAFFRON).

ADDIS depends on constants $w_{0},\lambda$ and $\tau$. $w_{0}$
represents the initial \`wealth’ of the procedure and satisfies
$0 \leq w_{0} \leq \alpha$. $\tau \in (0,1\rbrack$ represents the
threshold for a hypothesis to be selected for testing: p-values greater
than $\tau$ are implicitly ‘discarded’ by the procedure. Finally,
$\lambda \in \lbrack 0,\tau)$ sets the threshold for a p-value to be a
candidate for rejection: ADDIS will never reject a p-value larger than
$\lambda$.

The significance thresholds for ADDIS are chosen as follows:
$$\alpha_{t} = \min\{\lambda,{\widetilde{\alpha}}_{t}\}$$ where
$${\widetilde{\alpha}}_{t} = (\tau - \lambda)\lbrack w_{0}\gamma_{S^{t} - C_{0 +}} + \left( \alpha - w_{0} \right)\gamma_{S^{t} - \kappa_{1}^{*} - C_{1 +}} + \alpha\sum\limits_{j \geq 2}\gamma_{S^{t} - \kappa_{j}^{*} - C_{j +}}$$
and
$$\kappa_{j} = \min\{ i \in \lbrack t - 1\rbrack:\sum\limits_{k \leq i}1\{ p_{k} \leq \alpha_{k}\} \geq j\},\;\kappa_{j}^{*} = \sum\limits_{i \leq \kappa_{j}}1\{ p_{i} \leq \tau\},\; S^{t} = \sum\limits_{i < t}1\{ p_{i} \leq \tau\},\; C_{j +} = \sum\limits_{i = \kappa_{j} + 1}^{t - 1}1\{ p_{i} \leq \lambda\}$$

The default sequence of $\gamma$ for ADDIS is the same as for SAFFRON
given [here](#SAFFRON_gamma).

## FWER Control

### Alpha-spending

The Alpha-spending procedure controls the FWER for a potentially
infinite stream of p-values using a Bonferroni-like test. Given an
overall significance level $\alpha$, the significance thresholds are
chosen as $$\alpha_{i} = \alpha\gamma_{i}$$ where
$\sum_{i = 1}^{\infty}\gamma_{i} = 1$ and $\gamma_{i} \geq 0$. The
procedure strongly controls the FWER for arbitrarily dependent p-values.

Note that the procedure also controls the generalised familywise error
rate (k-FWER) for $k > 1$ if $\alpha$ is replaced by $\min(1,k\alpha)$.

The default sequence of $\gamma$ is the same as that for $\xi$ for LORD
given [here](#LORD_gamma).

### Online Fallback

The online fallback procedure of Tian & Ramdas (2019b) provides a
uniformly more powerful method than Alpha-spending, by saving the
significance level of a previous rejection. More specifically, online
fallback tests hypothesis $H_{i}$ at level
$$\alpha_{i} = \alpha\gamma_{i} + R_{i - 1}\alpha_{i - 1}$$ where
$R_{i} = 1\{ p_{i} \leq \alpha_{i}\}$ denotes a rejected hypothesis. The
procedure strongly controls the FWER for arbitrarily dependent p-values.

The default sequence of $\gamma$ is the same as that for $\xi$ for LORD
given [here](#LORD_gamma).

### ADDIS-spending

The ADDIS-spending procedure strongly controls the FWER for independent
p-values, and was proposed by Tian & Ramdas (2021). The procedure
compensates for the power loss of Alpha-spending, by including both
adapativity in the fraction of null hypotheses and the conservativeness
of nulls.

ADDIS depends on constants $\lambda$ and $\tau$, where $\lambda < \tau$.
Here $\tau \in (0,1)$ represents the threshold for a hypothesis to be
selected for testing: p-values greater than $\tau$ are implicitly
\`discarded’ by the procedure, while $\lambda \in (0,1)$ sets the
threshold for a p-value to be a candidate for rejection: ADDIS-spending
will never reject a p-value larger than $\lambda$.

Note that the procedure controls the generalised familywise error rate
(k-FWER) for $k > 1$ if $\alpha$ is replaced by $\min(1,k\alpha)$. Tian
and Ramdas (2019b) also presented a version for handling local
dependence, see the Section on Asynchronous testing below.

The default sequence of $\gamma$ for ADDIS-spending is the same as for
SAFFRON given [here](#SAFFRON_gamma).

## Accounting for dependent p-values

As noted above, the LORD, SAFFRON, ADDIS and ADDIS-spending procedures
assume independent p-values, while the LOND procedure is also valid
under positive dependencies (like the Benjamini-Hochberg method, see
below). Adjusted versions of LOND and LORD available for arbitrarily
dependent p-values. Alpha-spending and online fallback also control the
FWER and FDR for arbitrarily dependent p-values.

By way of comparison, the usual Benjamini-Hochberg method for
controlling the FDR assumes that the p-values are positively dependent
(PRDS). As an example, the PRDS is satisfied for multivariate normal
test statistics with a positive correlation matrix). See Benjamini &
Yekutieli (2001) for further technical details.
