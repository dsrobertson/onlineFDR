# The theory behind onlineFDR

## FDR Control

Consider a sequence of hypotheses $`H_1, H_2, H_3, \ldots`$ that arrive
sequentially in a stream, with corresponding $`p`$-values
$`(p_1, p_2, p_3, \ldots)`$. A testing procedure provides a sequence of
adjusted significance thresholds $`\alpha_i`$, with corresponding
decision rule:
``` math
R_i = \left\{\begin{array}{ccc}
1 & \text{if } p_i \leq \alpha_i & (\text{reject } H_i)\\
0 & \text{otherwise } & (\text{accept } H_i)
\end{array}\right\}
```

In *online* testing, the significance thresholds can only be functions
of the prior decisions,
i.e. $`\alpha_i = \alpha_i(R_1, R_2, \ldots, R_{i-1})`$.

Javanmard and Montanari (2015, 2018) proposed two procedures for online
control. The first is LOND, which stands for (significance) Levels based
On Number of Discoveries. The second is LORD, which stands for
(significance) Levels based On Recent Discovery. LORD was subsequently
extended by Ramdas *et al.* (2017). Ramdas *et al.* (2018) also proposed
the SAFFRON procedure, which provides an adaptive method of online FDR
control, which includes a variant of Alpha-investing. Finally, Tian and
Ramdas (2019) proposed the ADDIS procedure as an improvement of SAFFRON
in the presence of conservative nulls.

### LOND

The LOND procedure controls the FDR for independent $`p`$-values. Given
an overall significance level $`\alpha`$, we choose a sequence of
non-negative numbers $`\beta = (\beta_i)_{i \in \mathbb{N}}`$ such that
they sum to $`\alpha`$. The values of the adjusted significance
thresholds $`\alpha_i`$ are chosen as follows:
``` math
 \alpha_i = \beta_i (D(i-1) + 1) 
```
where $`D(n) = \sum_{i=1}^n R_i`$ denotes the number of discoveries
(i.e. rejections) in the first $`n`$ hypotheses tested.

A slightly modified version of LOND, with
$`\alpha_i = \beta_i \max(D(i-1), 1)`$ (Javanmard and Montanari, 2015,
Remark 2.3), also controls the FDR for positively dependent (PRDS)
$`p`$-values (Zrnic *et al.*, 2021). It is used by setting
`original = FALSE`. Their argument also shows that the original version
controls the FDR for PRDS $`p`$-values: if the $`i`$-th hypothesis is
rejected, the total number of discoveries is at least $`D(i-1) + 1`$,
the factor used in $`\alpha_i`$, and this factor can only increase when
$`p`$-values decrease, so Lemma 1 of Ramdas *et al.* (2019) applies.

LOND can be adjusted to also control FDR under arbitrarily dependent
$`p`$-values. To do so, it is modified with
$`\tilde{\beta}_i = \beta_i/H(i)`$ in place of $`\beta_i`$, where
$`H(i) = \sum_{j=1}^i \frac{1}{j}`$ is the $`i`$-th harmonic number.
Note that this leads to a substantial loss in power compared to the
unadjusted LOND procedure. The correction factor is similar to the
classical one used by Benjamini and Yekutieli (2001), except that in
this case the $`i`$-th hypothesis among $`N`$ is penalised by a factor
of $`H(i)`$ to give consistent results across time (as compared to a
factor $`H(N)`$ for the Benjamini and Yekutieli method).

The default sequence of $`\beta`$ is given by
``` math
\beta_j = C \alpha \frac{\log(\max(j, 2))}{j e^{\sqrt{\log j}}}
```
where $`C \approx 0.07908197`$ is chosen so that the $`\beta_j`$ sum to
$`\alpha`$. This is the sequence proposed by Javanmard and Montanari
(2018), equation 31, who give $`C \approx 0.07720838`$; with that value
the sum is about $`0.976\alpha`$.

### LORD

The LORD procedure controls the FDR for independent $`p`$-values. We
first fix a sequence of non-negative numbers
$`\gamma = (\gamma_i)_{i \in \mathbb{N}}`$ such that
$`\gamma_i \geq \gamma_j`$ for $`i \leq j`$ and
$`\sum_{i=1}^{\infty} \gamma_i = 1`$. At each time $`i`$, let $`\rho_i`$
be the last time a discovery was made before $`i`$:
``` math
\rho_i = \max \left\{ l \in \{1, \ldots, i-1\} : R_l = 1\right\}
```
with $`\rho_i = 0`$ if no discovery has been made before $`i`$. Also let
$`\kappa_j`$ be the time of the $`j`$-th discovery:
``` math
\kappa_j = \min \left\{ i \in \mathbb{N} : \sum_{l \leq i} R_l \geq j \right\}
```

LORD depends on constants $`w_0`$ and $`b_0`$, where $`w_0 \geq 0`$
represents the initial ‘wealth’ of the procedure and $`b_0 > 0`$
represents the ‘payout’ for rejecting a hypothesis. We require
$`w_0+b_0 \leq \alpha`$ for FDR control to hold.

Javanmard and Montanari (2018) presented three different versions of
LORD, which have different definitions of the adjusted significance
thresholds $`\alpha_i`$. Versions 1 and 2 have since been superseded by
the LORD++ procedure of Ramdas *et al.* (2017), so we do not describe
them here.

- **LORD++**: The significance thresholds for LORD++ are chosen as
  follows:
  ``` math
  \alpha_i = \gamma_i w_0 + (\alpha - w_0) \gamma_{i-\kappa_1} +
  \alpha \sum_{j : \kappa_j < i, \kappa_j \neq \kappa_1} \gamma_{i - \kappa_j}
  ```

- **LORD 3**: The significance thresholds depend on the time of the last
  discovery and the wealth accumulated at that time, with
  ``` math
  \alpha_i  = \gamma_{i - \rho_i} W(\rho_i)
  ```
  where $`\rho_i`$ is as defined above (so $`\rho_1 = 0`$). Here
  $`\{W(j)\}_{j \geq 0}`$ represents the ‘wealth’ available at time
  $`j`$, and is defined recursively:

``` math
\begin{aligned}
W(0) &= w_0 \\
W(j) &= W(j-1) - \alpha_j + b_0 R_j
\end{aligned}
```

- **D-LORD**: This is equivalent to the LORD++ procedure with
  discarding. Given a discarding threshold $`\tau \in (0,1]`$ and
  initial wealth $`w_0 \leq \tau\alpha`$ the significance thresholds are
  chosen as follows:
  ``` math
  \alpha_t = \min\{\tau, \tilde{\alpha}_t\}
  ```
  where
  ``` math
  \tilde{\alpha}_t = w_0 \gamma_{S^t} +
  (\tau\alpha - w_0)\gamma_{S^t - \kappa_1^*} +
  \tau\alpha \sum_{j \geq 2} \gamma_{S^t - \kappa_j^*}
  ```
  and
  ``` math
  \kappa_j = \min\{i \in [t-1] : \sum_{k \leq i}
  1 \{p_k \leq \alpha_k\} \geq j\}, \;
  \kappa_j^* = \sum_{i \leq \kappa_j} 1 \{p_i \leq \tau \}, \;
  S^t = \sum_{i < t} 1 \{p_i \leq \tau \}
  ```
  Here the sequence $`\gamma`$ is indexed from 0, as in Tian and Ramdas
  (2019), so that the first element of `gammai` is $`\gamma_0`$.

LORD++ is an instance of a monotone rule, and provably controls the FDR
for independent p-values provided $`w_0 \leq \alpha`$. LORD 3 is a
non-monotone rule, and FDR control is only demonstrated empirically. In
some scenarios with large $`N`$, LORD 3 will have a slightly higher
power than LORD++ (see Robertson and Wason, 2018), but since it is a
non-monotone rule we would recommend using LORD++ (which is the
default), especially since it also has a provable guarantee of FDR
control.

In all versions except ‘dep’, the default sequence of $`\gamma`$ is
given by
``` math
\gamma_j = C \frac{\log(\max(j, 2))}{j e^{\sqrt{\log j}}}
```
where $`C \approx 0.07908197`$ is chosen so that the $`\gamma_j`$ sum
to 1. This is the sequence proposed by Javanmard and Montanari (2018),
equation 31, who give $`C \approx 0.07720838`$; with that value the sum
is about 0.976.

Javanmard and Montanari (2018) also proposed an adjusted version of LORD
that is valid for arbitrarily *dependent* p-values. Similarly to LORD 3,
the adjusted significance thresholds are set equal to
``` math
 \alpha_i = \xi_i W(\rho_i)
```
where (assuming $`w_0 \leq b_0`$),
$`\sum_{j=1}^{\infty} \xi_j (1 + \log(j)) \leq \alpha / b_0`$ (Javanmard
and Montanari, 2018, example 3.7). If $`w_0 > b_0`$, the condition is
instead $`\sum_{j=1}^{\infty} (w_0 + b_0 \log(j)) \xi_j \leq \alpha`$,
which follows from their Theorem 3.6; the default sequence is then
$`\xi_j = C'/(j \log(\max(j, 2))^3)`$, with $`C'`$ chosen so that this
sum equals $`\alpha`$.

The default sequence of $`\xi`$ is given by
``` math
 \xi_j  = \frac{C}{j \log(\max(j, 2))^3}
```
where $`C \approx 0.1393066`$ is chosen so that
$`\sum_{j=1}^{\infty} \xi_j (1 + \log(j)) = 1 \leq \alpha / b_0`$. (The
deprecated function `LORDdep` multiplies this sequence by
$`\alpha / b_0`$.)

Note that allowing for dependent p-values can lead to a substantial loss
in power compared with the LORD procedures described above.

### SAFFRON

The SAFFRON procedure controls the FDR for independent p-values, and was
proposed by Ramdas *et al.* (2018). The algorithm is based on an
estimate of the proportion of true null hypotheses. More precisely,
SAFFRON sets the adjusted test levels based on an estimate of the amount
of alpha-wealth that is allocated to testing the true null hypotheses.

SAFFRON depends on constants $`w_0`$ and $`\lambda`$, where $`w_0`$
satisfies $`0 \leq w_0 \leq \alpha`$ and represents the initial ‘wealth’
of the procedure, and $`\lambda \in (0,1)`$ represents the threshold for
a ‘candidate’ hypothesis. A ‘candidate’ refers to p-values no larger
than $`\lambda`$, since SAFFRON will never reject a p-value larger than
$`\lambda`$. These candidates can be thought of as the hypotheses that
are more likely to be non-null.

The SAFFRON procedure runs as follows:

1.  At each time $`t`$, define the number of candidates after the
    $`j`$-th rejection as
    ``` math
     C_{j+} = C_{j+}(t) = \sum_{i = \kappa_j + 1}^{t-1} C_i
    ```
    where $`C_t = 1\{p_t \leq \lambda \}`$ is the indicator for
    candidacy, $`\kappa_j`$ is the time of the $`j`$-th rejection (as
    defined for LORD) and $`\kappa_0 = 0`$.

2.  SAFFRON starts with
    $`\alpha_1 = \min\{(1 - \lambda)\gamma_1 w_0, \lambda\}`$.
    Subsequent test levels are chosen as
    $`\alpha_t = \min\{ \lambda, \tilde{\alpha}_t\}`$, where
    ``` math
    \tilde{\alpha}_t = (1 - \lambda) [w_0 \gamma_{t-C_{0+}} +
    (\alpha - w_0)\gamma_{t-\kappa_1-C_{1+}} +
    \alpha \sum_{j \geq 2}
    \gamma_{t - \kappa_j- C_{j+}}]
    ```

The default sequence of $`\gamma`$ for SAFFRON is given by
$`\gamma_j \propto j^{-1.6}`$.

### Alpha-investing

Ramdas *et al.* (2018) proposed a variant of the Alpha-investing
algorithm of Foster and Stine (2008) that guarantees FDR control for
independent p-values (but not for positively dependent p-values in
general; see [Accounting for dependent
p-values](#accounting-for-dependent-p-values)). This procedure uses
SAFFRON’s update rule with the constant $`\lambda`$ replaced by a
sequence $`\lambda_i = \alpha_i`$. This is also equivalent to using the
ADDIS algorithm (see below) with $`\tau = 1`$ and
$`\lambda_i = \alpha_i`$.

### ADDIS

The ADDIS procedure controls the FDR for independent p-values, and was
proposed by Tian and Ramdas (2019). The algorithm compensates for the
power loss of SAFFRON with conservative nulls, by including both
adaptivity in the fraction of null hypotheses (like SAFFRON) and the
conservativeness of nulls (unlike SAFFRON).

ADDIS depends on constants $`w_0, \lambda`$ and $`\tau`$. $`w_0`$
represents the initial ‘wealth’ of the procedure and satisfies
$`0 \leq w_0 \leq \alpha`$. $`\tau \in (0,1]`$ represents the threshold
for a hypothesis to be selected for testing: p-values greater than
$`\tau`$ are implicitly ‘discarded’ by the procedure. Finally,
$`\lambda \in (0,\tau)`$ sets the threshold for a p-value to be a
candidate for rejection: ADDIS will never reject a p-value larger than
$`\lambda`$.

The significance thresholds for ADDIS are chosen as follows:
``` math
\alpha_t = \min\{\lambda, \tilde{\alpha}_t\}
```
where
``` math
\tilde{\alpha}_t = (\tau - \lambda)[w_0 \gamma_{S^t-C_{0+}} + 
(\alpha - w_0)\gamma_{S^t - \kappa_1^*-C_{1+}} +
\alpha \sum_{j \geq 2} \gamma_{S^t - \kappa_j^* - C_{j+}}]
```
and
``` math
\kappa_j = \min\{i \in [t-1] : \sum_{k \leq i}
1 \{p_k \leq \alpha_k\} \geq j\}, \;
\kappa_j^* = \sum_{i \leq \kappa_j} 1 \{p_i \leq \tau \}, \;
S^t = \sum_{i < t} 1 \{p_i \leq \tau \}, \;
C_{j+} = \sum_{i = \kappa_j + 1}^{t-1} 1\{p_i \leq \lambda\}
```
Here the sequence $`\gamma`$ is indexed from 0, as in Tian and Ramdas
(2019), so that the first element of `gammai` is $`\gamma_0`$.

The default sequence of $`\gamma`$ for ADDIS is the same as for SAFFRON
given [here](#SAFFRON_gamma).

## FWER Control

### Alpha-spending

The Alpha-spending procedure controls the FWER for a potentially
infinite stream of p-values using a Bonferroni-like test. Given an
overall significance level $`\alpha`$, the significance thresholds are
chosen as
``` math
\alpha_i = \alpha \gamma_i
```
where $`\sum_{i=1}^{\infty} \gamma_i = 1`$ and $`\gamma_i \geq 0`$. The
procedure strongly controls the FWER for arbitrarily dependent p-values.

Note that the procedure also controls the generalised familywise error
rate (k-FWER) for $`k > 1`$ if $`\alpha`$ is replaced by
$`\min(1,k\alpha)`$.

The default sequence of $`\gamma`$ is the same as that for $`\gamma`$
for LORD given [here](#LORD_gamma).

### Online Fallback

The online fallback procedure of Tian and Ramdas (2021) provides a
uniformly more powerful method than Alpha-spending, by saving the
significance level of a previous rejection. More specifically, online
fallback tests hypothesis $`H_i`$ at level
``` math
\alpha_i = \alpha \gamma_i + R_{i-1} \alpha_{i-1}
```
for $`i \geq 2`$, with $`\alpha_1 = \alpha\gamma_1`$, where
$`R_i = 1\{p_i
\leq \alpha_i\}`$ denotes a rejected hypothesis. The procedure strongly
controls the FWER for arbitrarily dependent p-values.

The default sequence of $`\gamma`$ is the same as that for $`\gamma`$
for LORD given [here](#LORD_gamma).

### ADDIS-spending

The ADDIS-spending procedure strongly controls the FWER for independent
p-values, and was proposed by Tian and Ramdas (2021). The procedure
compensates for the power loss of Alpha-spending, by including both
adaptivity in the fraction of null hypotheses and the conservativeness
of nulls.

ADDIS-spending depends on constants $`\lambda`$ and $`\tau`$, where
$`\lambda < \tau`$. Here $`\tau \in (0,1]`$ represents the threshold for
a hypothesis to be selected for testing: p-values greater than $`\tau`$
are implicitly ‘discarded’ by the procedure, while $`\lambda \in (0,1)`$
sets the threshold for a p-value to be a candidate for rejection:
ADDIS-spending will never reject a p-value larger than $`\lambda`$.

Note that the procedure controls the generalised familywise error rate
(k-FWER) for $`k > 1`$ if $`\alpha`$ is replaced by $`\min(1,k\alpha)`$.
Tian and Ramdas (2021) also presented a version for handling local
dependence (see
[`?ADDIS_spending`](https://dsrobertson.github.io/onlineFDR/reference/ADDIS_spending.md)).

The default sequence of $`\gamma`$ for ADDIS-spending is the same as for
SAFFRON given [here](#SAFFRON_gamma).

## Accounting for dependent p-values

The LORD, SAFFRON, Alpha-investing, ADDIS, ADDIS-spending and
ADDIS-exhaustive procedures assume that the null p-values are
independent of each other and of the non-null p-values. Unlike for the
Benjamini-Hochberg procedure (see below), their guarantees do not extend
to positively dependent (PRDS) p-values in general, even when the
dependence is confined to blocks of consecutive tests. Fisher (2024)
shows that LORD and SAFFRON still control the FDR under positive
dependence if each null p-value is conditionally super-uniform given the
information used to set its test level. For well-powered tests this
essentially requires each p-value to be independent of that information
(local dependence in Fisher’s sense), for example when the levels for a
block of tests are set only from the results of earlier, independent
blocks.

Without such a condition the FDR can exceed the nominal level. For
example, with one-sided tests of normal means whose test statistics
share a common factor, either across the whole sequence or within blocks
of tests ($`\alpha = 0.05`$):

- SAFFRON with the default settings, all hypotheses null: the FDR is
  0.051 for 30 tests with pairwise correlation 0.5, and 0.068 for 100
  tests and 0.075 for 1000 tests with pairwise correlation 0.7. The
  inflation is larger at smaller $`\alpha`$: for $`\alpha = 0.01`$ and
  1000 tests the FDR is about $`4\alpha`$.
- ADDIS with the default settings: the FDR is about 0.085 when 20 strong
  signals are followed by 80 nulls and all test statistics have pairwise
  correlation 0.9.
- Alpha-investing with the default settings: in the same setting the FDR
  is about 0.15, three times the nominal level.
- LORD++ with the default settings: the FDR is about 0.063 for 10,500
  tests in 50 independent blocks, each of 10 strong signals (p-values of
  essentially 0) followed by 200 nulls whose test statistics have
  pairwise correlation 0.99 within the block. With $`\gamma_j = 2^{-j}`$
  and blocks of 10 strong signals followed by 50 such nulls (3,000
  tests), the FDR is about 0.135.

For SAFFRON, ADDIS and Alpha-investing the reason is that alpha-wealth
is not spent on some p-values: those below $`\lambda`$ for SAFFRON,
those outside $`(\lambda, \tau]`$ for ADDIS, and the rejected ones for
Alpha-investing. When a shared effect makes all the p-values small, the
test levels stop decreasing exactly when the null p-values are also
small. For LORD++, each rejection raises the levels of the tests that
follow it, so when a block of correlated null p-values is small, one
false rejection makes the next ones more likely and the rejections
cascade.

For p-values that may be positively dependent across the sequence, LOND
controls the FDR under the PRDS condition (see above), and so does
BatchPRDS when the p-values are independent across batches. Adjusted
versions of LOND (`dep = TRUE`) and LORD (`version = 'dep'`) control the
FDR for arbitrarily dependent p-values, while LONDstar, LORDstar and
SAFFRONstar (`version = 'dep'`) control the modified FDR (mFDR) under
local dependence. Alpha-spending and online fallback control the FWER
(and hence the FDR) for arbitrarily dependent p-values, and
ADDIS-spending (`dep = TRUE`) controls the FWER under local dependence.

By way of comparison, the usual Benjamini-Hochberg method for
controlling the FDR assumes that the p-values are positively dependent
(PRDS). As an example, the PRDS is satisfied for one-sided p-values from
multivariate normal test statistics whose correlations are all
non-negative. See Benjamini and Yekutieli (2001) for further technical
details.

## References

Benjamini, Y. and Yekutieli, D. (2001). The control of the false
discovery rate in multiple testing under dependency. *The Annals of
Statistics*, 29(4):1165-1188.

Fisher, A. (2024). Online false discovery rate control for LORD++ and
SAFFRON under positive, local dependence. *Biometrical Journal*,
66(1):2300177.

Foster, D. and Stine, R. (2008). $`\alpha`$-investing: a procedure for
sequential control of expected false discoveries. *Journal of the Royal
Statistical Society (Series B)*, 70(2):429-444.

Javanmard, A. and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, <https://arxiv.org/abs/1502.06197>.

Javanmard, A. and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Ramdas, A., Barber, R.F., Wainwright, M.J. and Jordan, M.I. (2019). A
unified treatment of multiple testing with prior knowledge using the
p-filter. *The Annals of Statistics*, 47(5):2790-2821.

Ramdas, A., Yang, F., Wainwright, M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Ramdas, A., Zrnic, T., Wainwright, M.J. and Jordan, M.I. (2018).
SAFFRON: an adaptive algorithm for online control of the false discovery
rate. *Proceedings of the 35th International Conference on Machine
Learning*, 80:4286-4294.

Robertson, D.S. and Wason, J.M.S. (2018). Online control of the false
discovery rate in biomedical research. *arXiv preprint*,
<https://arxiv.org/abs/1809.07292>.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems 32*, 9388-9396.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods in Medical Research*, 30(4):976-993.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*,
22(33):1-39.
