---
bibliography: [references.bib]
---

# Introduction

Species interactions drive many processes in evolution and ecology. A better
understanding of species interactions is an imperative to understand the
evolution of life on Earth, to mitigate the impacts of anthropogenic change on
biodiversity [@Makiola2020KeyQue], and for predicting zoonotic spillover of
disease to prevent future pandemics [@Becker2021OptPre]. At the moment we lack
sufficient data to meet these challenges [@Poisot2021GloKno], largely because
species interactions are hard to sample [@Jordano2016SamNet]. Over the past few
decades biodiversity data has become increasingly available through remotely
collected data and adoption of open data practices
[@Stephenson2020TecAdv;@Kenall2014OpeFut]. Still, interaction data remains
relatively scarce because sampling typically requires human observation. This
induces a constraint on the amount, spatial scale, and temporal frequency of
resulting data that it is feasible to collect by humans. Many crowdsourced
methods for biodiversity data aggregation (e.g. GBIF, eBird) still rely on
automated identification of species, which does not easily generalize to
interaction sampling. There is interest in using remote methods for interaction
sampling, which primarily detect co-occurrence and derive properties like
species avoidance from this data [@Niedballa2019AssAna]. However, co-occurrence
itself is not necessarily indicative of an interaction [@Blanchet2020CooNot].
This is an example of semantic confusion around the word “interaction”---for
example one might consider competition a type of species interaction, even
though it is marked by a lack of co-occurrence between species, unlike other
types of interactions, like predation or parasitism, which require both species
to be together at the same place and time. Here we consider interaction in the
latter sense, where two species have fitness consequences on one-another if (and
only if) they are in the sample place at the same time. In addition, here we
only consider direct (not higher-order) interactions. 

We cannot feasibly observe all (or even most) of the interactions that occur in
an ecosystem. This means we can be confident two species actually interact if we
have a record of it (assuming they are correctly identified), but not at all
confident that a pair of species _do not_ interact if we have _no record_ of
those species observed together. In other words, it is difficult to distinguish
_true-negatives_ (two species never interact) from _false-negatives_ (two
species interact sometimes, but we do not have a record of this interaction). For a concrete
example of a false-negative in a food web, see @fig:concept. Because even the
most highly sampled systems will still contain false-negatives, there is
increasing interest in combining species-level data (e.g. traits, abundance,
range, phylogenetic relatedness, etc.) to build models to predict interactions
between species we haven't observed together before [@Strydom2021RoaPre].
However, the noise of false-negatives could impact the efficacy of our
predictive models and have practical consequences for answering questions about
interactions [@deAguiar2019RevBia]. This data constraint is amplified as the
interaction data we have is geographically biased toward the usual suspects
[@Poisot2021GloKno]. We therefore need a statistical approach to assessing these
biases in the observation process and their consequences for our understanding
of interaction networks. 

The importance of _sampling effort_ and its impact on resulting ecological data
has produced a rich body of literature.The recorded number of species in a
dataset or sample depends on the total number of observations
 [@Willott2001SpeAcc; @Walther1995SamEff], as do estimates of population
abundance [@Griffiths1998SamEff]. This relationship between sampling effort and
spatial coverage and species detectability has motivated more quantitatively
robust approaches to account for error in sampling data in many contexts: to
determine if a given species is extinct [@Boakes2015InfSpe], to determine
sampling design [@Moore2016OptEco], and to measure species richness across large
scales [@Carlson2020WhaWou]. In the context of interactions, an initial concern
was the compounding effects of limited sampling effort combined with the
amalgamation of data (across both study sites, time of year, and taxonomic
scales) could lead any empirical set of observations to inadequately reflect the
reality of how species interact [@Paine1988RoaMap] or the structure of the
network as a whole [@Martinez1999EffSam; @McLeod2021SamAsy]. @Martinez1999EffSam
showed that in a plant-endophyte trophic network, network connectance is robust
to sampling effort, but this was done in the context of a system for which
observation of 62,000 total interactions derived from 164,000 plant-stems was
feasible. In some systems (e.g. megafauna food-webs) this many observations is
either impractical or infeasible due to the absolute abundance of the species in
question. 


The intrinsic properties of ecological communities create several challenges for
sampling: first, species are not observed with equal probability---we are much
more likely to observe a species of high abundance than one of very low
abundance [@Poisot2015SpeWhy]. @Canard2012EmeStr presents a null model of
food-web structure where species encounter one-another in proportion to each
species’ relative-abundance. This assumes that there are no associations in
species co-occurrence due to an interaction (perhaps because this interaction is
"important" for both species; @Cazelles2016TheSpe), but in this paper we later
show increasing strength of associations leads to increasing probability of
false-negatives in interaction data, and that these positive associations are
common in existing network data. Second, observed co-occurrence is often
equated with meaningful interaction strength, but this is not necessarily the
case [@Blanchet2020CooNot]---a true “non-interaction" would require that neither
of two species, regardless of whether they co-occur, ever exhibit any meaningful
effect on the fitness of the other. So, although co-occurrence is not directly
indicative of an interaction, it _is_ a precondition for an interaction. 

Here, we illustrate how our confidence that a pair of species never interacts
highly depends on sampling effort. We demonstrate how the realized
false-negative-rate of interactions is related to the relative abundance of the
species pool, and introduce a method to produce a null estimate of the
false-negative-rate given total sampling effort (the total count of all
interactions seen among all species-pairs) and a method for including
uncertainty into model predictions of interaction probabilities to account for
observation error. We then confront these models with data, by showing that
positive associations in co-occurrence data can increase the realized number of
false-negatives and by showing these positive associations are rampant in
network datasets. We conclude by recommending that the simulation of sampling
effort and species occurrence can and should be used to help design surveys of
species interaction diversity [@Moore2016OptEco], and by advocating use of null
models like those presented here as a tool for both guiding design of surveys of
species interactions and for including detection error into predictive models.


![This conceptual example considers a sample of the trophic community of bears, wolves, salmon (pink fish), pike (yellow fish), berry trees, and aspen trees. The true metaweb (all realized interactions across the entire spatial extent) is shown on the left. In the center is what a hypothetical ecologist samples at each site. Notice that although bears are observed co-occurring with both salmon and pike, there was never a direct observation of bears eating pike, even though they actually do. Therefore, this interaction between bears and pike is a false-negative.](./figures/concept.png){#fig:concept}

# Accounting for false-negatives in species interactions

In this section, we demonstate how difference in relative-abundance can lead to
many false-negatives in interaction data. We also introduce a method for
producing a null estimate of the false-negative-rate in datasets via simulation,
and a method for incorporating uncertainty directly into predictions of species
interactions to account for observation error.

## How many observations of a non-interaction do we need to be confident it's a true negative?

We start with a naive model of interaction detection: we assume that every
interacting pair of species is incorrectly observed as not-interacting with an
independent and fixed probability, which we denote $p_{fn}$ and subsequently refer
to as the False-Negative-Rate (FNR). If we observe the same species
not-interacting $N$ times, then the probability of a true-negative (denoted
$p_{tn}$) is given by $p_{tn}=1-(p_{fn})^N$. This relation (the
probability-mass-function of geometric distribution, a special case of the
negative-binomial distribution) is shown in @fig:geometric(A) for varying values
of $p_{fn}$ and illustrates a fundamental link between our ability to reliably
say an interaction doesn't exist---$p_{tn}$---and the number of times $N$ we
have observed a given species. In addition, note that there is no non-zero
$p_{fn}$ for which we can ever _prove_ that an interaction does not exist---no
matter how many observations of non-interactions $N$ we have, $p_{tn}<1$.

From @fig:geometric(A) it is clear that the more often we see two species
co-occurring, but _not interacting_, the more likely the interaction is a
true-negative. This has several practical consequences: first it means negatives
taken outside the overlap of the range of each species aren’t informative
because co-occurrence was not possible, and therefore neither was an
interaction. Second, we can use this relation to compute the expected number of
total observations needed to obtain a "goal" number of observations of a
particular pair of species (@fig:geometric(B)). As an example, if we hypothesize
that $A$ and $B$ do not interact, and we want to see species $A$ and $B$ both
co-occurring and _not interacting_ 10 times to be confident this is a true
negative, then we need an expected 1000 observations of all species if the
relative abundances of $A$ and $B$ are both $0.1$.

Because the true FNR is latent, we can never actually be sure what the actual
number of false-negatives in our data---however, we can use simulation to
estimate the FNR for datasets of a given size using neutral models of
observation. If some of the “worst-case” FNRs presented in @fig:geometric(A)
seem unrealistically high, considering that species are observed in proportion
to their relative abundance. In the next section we demonstrate that the
distribution of abundance in ecosystems can lead to very high realized values of
FNR ($p_{fn}$) simply as an artifact of sampling effort.


![**(A)** The probability that an observed interaction is a true negative
(y-axis) given how many times it has been sampled as a non-interaction (x-axis).
Each color reflects a different value of $p_{fn}$, the false-negative-rate
(FNR)---this is effectively the cumulative distribution function (cdf) of the
geometric distribution. **(B)** The expected number of total observations needed
(colors) to observe 10
co-occurrences between a species with relative abundance $P(A)$ (x-axis) and a
second species with relative abundance $P(Y)$. **(C)**: false-negative-rate
(y-axis) as a function of total sampling effort (x-axis) and network size,
computed using the method described above. For 500 independent draws from the
niche model (@Williams2000SimRul) at varying levels of species richness (colors)
with connectance drawn according to the flexible-links model
(@MacDonald2020RevLin) as described in the main text. For each draw from the
niche model, 200 sets of 1500 observations are simulated, for which the
mean false-negative-rate at each observation-step is computed. Means denoted
with points, with 1 in the first shade and 2 in the second. **(D)**: Same as
**(C)**, except using empirical food webs from Mangal database, where richness.
The outlier on **(D)** is a 714 species
food-web.](./figures/fig1.png){#fig:geometric}


## False-negatives as a product of relative abundance

We now show that the realized FNR changes drastically with sampling effort due
to the intrinsic variation of the abundance of individuals of each species
within a community. We do this by simulating the process of observation of
species interactions, applied both to 243 empirical food webs from the Mangal
database [@Banville2021ManJl] and random food-webs generated using the niche
model, a simple generative model of food-web structure that accounts for
allometric scaling [@Williams2000SimRul]. Our neutral model of observation
assumes each observed species is drawn in proportion to each species' abundance
at that place and time. The abundance distribution of a community can be
reasonably-well described by a log-normal distribution [@Volkov2003NeuThe]. In
addition to the log-normal distribution, we also tested the case where the
abundance distribution is derived from power-law scaling $Z^{(log(T_i)-1)}$
where $T_i$ is the trophic level of species $i$ and $Z$ is a scaling coefficient
[@Savage2004EffBod], which yields the same qualitative behavior. The practical
consequence of abundance distributions spanning many orders of magnitude of
abundance is that observing two "rare" species interacting requires two low
probability events: observing two rare species _at the same time_.

To simulate the process of observation, for an ecological network $M$ with $S$
species, we sample abundances for each species from a standard-log-normal
distribution. For each true interaction in the adjacency matrix $M$ (i.e.
$M_{ij}=1$) we estimate the probability of observing both species $i$ and $j$ at
a given place and time by simulating $n$ observations of all individuals of any
a species, where the species of the individual observed at the
$\{1,2,\dots,n\}$-th observation is drawn from the generated log-normal
distribution of abundances. For each pair of species $(i,j)$, if both $i$ and
$j$ are observed within the n-observations, the interaction is tallied as a true
positive if $M_{ij}=1$. If only one of $i$ or $j$ are observed---but not
both---in these $n$ observations, but $M_{ij}=1$, this is counted as a
false-negative, and a true-negative otherwise. For each pair of species $(i,j)$,
if both $i$ and $j$ are observed within the n-observations, the interaction is
tallied as a true positive if $M_{ij}=1$. If only one of $i$ or $j$ are
observed---but not both---in these $n$ observations, but $M_{ij}=1$, this is
counted as a false-negative, and a true-negative otherwise ($M_{ij} = 0$). This
process is illustrated conceptually in @fig:resampling_concept(A). 


In @fig:geometric(C) we see this model of observation applied to niche model
networks across varying levels of species richness, and in @fig:geometric(D) the
observation model applied to Mangal food webs. For all niche model simulations
in this manuscript, for a given number of species $S$ the number of interactions
is drawn from the flexible-links model fit to Mangal data
[@MacDonald2020RevLin], effectively drawing the number of interactions $L$ for a
random niche model food-web as

$$L \sim  \text{BetaBinomial}(S^2-S+1,\mu\phi, 1-\mu\phi)$$

where the maximum _a posteriori_ (MAP) estimate of $(\mu, \phi)$ applied to Mangal data from
[@MacDonald2020RevLin] is $(\mu=0.086, \phi=24.3)$. All simulations were done
with 500 independent replicates of unique niche model networks per unique number
of observations $n$. All analyses presented here are done in Julia v1.8
[@Bezanson2015JulFre] using both EcologicalNetworks.jl v0.5 and Mangal.jl v0.4
[@Banville2021ManJl] and are hosted on \href{https://github.com/gottacatchenall/ms_false_negatives/tree/main/src}{Github}). Note that the
empirical data, for the reasons described above, very likely already contains
many false-negatives, we'll revisit this issue in the final section.

From @fig:geometric(C) it is evident that the number of species considered in a
study is inseparable from the false-negative-rate in that study, and this effect
should be taken into account when designing samples of ecological networks in
the future. We see a similar qualitative pattern in  @fig:geometric(D) where the
FNR drops off quickly as a function of observation effort, mediated by total
richness. The practical consequence of the bottom row of @fig:geometric is
whether the total number of observations of all species (the x-axis) for the
threshold FNR we deem acceptable (the y-axis) is feasible. This raises two
points: first, empirical data on interactions are subject to the practical
limitations of funding and human-work hours, and therefore existing data tend to
fall on the order of hundreds or thousands observations of individuals per site.
Clear aggregation of data on sampling effort has proven difficult to find and a
meta-analysis of network data and sampling effort seems both pertinent and
necessary, in addition to the effects of aggregation of interactions across
taxonomic scales [@Gauzens2013FooAgg; @Giacomuzzo2021FooWeb]. This inherent
limitation on in-situ sampling means we should optimize where we sample across
space so that for a given number of samples, we obtain the maximum information
possible.  Second, what is meant by “acceptable” FNR? This raises the question:
does a shifting FNR lead to rapid transitions in our ability inference and
predictions about the structure and dynamics of networks, or does it produce a
roughly linear decay in model efficacy? We explore this in the next section.

We conclude this section by advocating for the use of neutral models similar to
above to generate expectations about the number of false-negatives in a data set
of a given size. This could prove fruitful both for designing surveys of
interactions but also because we may want to incorporate models of imperfect
detection error into predictive interactions models, as @Joseph2020NeuHie does
for species occurrence modeling. Additionally, we emphasize that one must
consider the context for sampling---is the goal to detect a particular species
(as in @fig:geometric(C)), or to get a representative sample of interactions
across the species pool? These arguments are well-considered when sampling
individual species [@Willott2001SpeAcc], but have not yet been adopted for
designing samples of communities.

## Including observation error in interaction predictions 

Here we show how to incorporate uncertainty into model predictions of
interaction probability to account for imperfect observation (both
false-negatives and false-positives). Models for interaction prediction
typically yield a probability of interaction between each pair of species,
$p_{ij}$. When these are considered with uncertainty, it is usually
model-uncertainty, e.g. the variance in the interaction probability prediction
across several cross-validation folds, where the data is split into training and
test sets several times.  The method we introduce adjusts the value of a model's
predictions to produce a distribution of interaction probabilities, which are
adjusted by a given false-negative-rate $p_{fn}$ and false-positive-rate
$p_{fp}$ (outlined in figure @fig:resampling_concept). We describe first how to
sample from this distribution of adjusted interaction probabilities via
simulation, and show that this distribution can be well-approximated
analytically. 

![(A) The process for estimating the false-negative-rate (FNR) for an
interaction dataset consisting of $N$ total observed interactions. (B) The
method for resampling interaction probability based on estimates of
false-negative and false-positive rates. (C) The method for interaction
probability resampling applied to three mammals and three parasites from the
@Hadfield2014TalTwo dataset. The original probability $p_{ij}$ is indicated with
a vertical dashed line. The histogram is simulated from the resampling process,
and the line indicates the gaussian approximation to this distribution.](./figures/uncertainty_sampler.png){#fig:resampling_concept}

We then consider the output prediction from an arbitrary prediction model, which
is the probability $p_{ij}$ that two species $i$ and $j$ interact. To get an
estimate of $p_{ij}$ that accounts for observation error, we resample the
probability of each interaction $p_{ij}$ by simulating a set of several
'particles', where each particle is a realization of an interaction occurring
(either true or false with probabilities $p_{ij}$ and $1-p_{ij}$ respectively)
and then being correctly observed with probabilities given by $p_{fp}$ and
$p_{fn}$ to yield a single boolean outcome for each particle (“Resampling”
within @fig:resampling_concept (B)). Across of many particles, the
resulting frequency of ‘true’ outcomes is a single resample of the interaction
probability $p_{ij}^*$. Across several samples each of several particles, this
forms a distribution of probabilities which are adjusted by the true and false
negative rates.

There is also an analytic way to approximate this distribution using the normal
approximation to binomial. As a reminder, as the total number of samples $N$
from a binomial distribution for $n$ trials with success probability $p$ from approaches
infinity, the sum of total successes across all samples approaches a normal
distribution with mean $np$ and variance $np(1-p)$. We can use this to correct
the estimate $p_{ij}$ based on the expected false-negative-rate $p_{fn}$ and
false-positive rate $p_{fp}$ to obtain the limiting distribution as the number
of resamples approaches infinity for the resampled $p_{ij}^*$ for a given number
of particles $n_p$. We do this by first adjusting for the rates of observation
error to get the mean resampled probability, $\mathbb{E}[{p_{ij}^*}]$, as

$$
\mathbb{E}[{p_{ij}^*}] = p_{ij}(1-p_{fp})+ (1-p_{ij})p_{fn}
$$

which is obatinable by definition (supp 1.)

For notation, here we refer to a normal distribution with mean $\mu$ and standard-deviation
$\sigma$ as $\mathcal{N}(\mu,\sigma)$

Then yields the normal approximation

$$
\sum_{i=1}^{n_p} p_{ij}^* \sim \mathcal{N}\bigg(n_p \cdot \mathbb{E}[p_{ij}^*], \sqrt{n_p\mathbb{E}[p_{ij}^*] (1- \mathbb{E}[p_{ij}^*])}\bigg)
$$

which then can be converted back to a distribution of frequency of successes to
yield the final approximation

$$p_{ij}^* \sim \mathcal{N}\bigg( \mathbb{E}[p_{ij}^*] , \sqrt{\frac{\mathbb{E}[p_{ij}^*]
(1-\mathbb{E}[p_{ij}^*] )}{n_p}} \bigg)$$ 
{#eq:eq1}

We can then further truncate this distribution to remain on the interval
$(0,1)$, as the output is a probability, although in practice often the
probability mass outside $(0,1)$ is extremely low except for $p_{ij}$ values
very close to 0 or 1.  As an example case study, we use a
boosted-regression-tree to predict interactions in a host-parasite network
[@Hadfield2014TalTwo] (with features derived in the same manner as
@Strydom2021RoaPre derives features on this data) to produce a set of
interaction predictions. We then applied this method to a set of a few resampled
interaction probabilities between mammals and parasite species shown in figure
@fig:resampling_concept(C). 

Why is this useful? For one, this analytic method avoids the extra computation
required by simulating samples from this distribution directly. Further, it
enables the extension of the natural analogue between $n_p$ (the number of
particles) and the number of observations of co-occurrence for a given pair of
species---the fewer the particles, the higher the variance of the resulting
approximation. The normal approximation is undefined for 0 particles (i.e. 0
observations co-occurrence), although as $n_p$ approaches 0 the approximated
normal (once truncated) approaches the uniform distribution on the interval
$(0,1)$, the maximum entropy distribution where we have no information about the
possibility of an interaction. 

This also has implications for what we mean by ‘uncertainty’ in interaction
predictions. A model’s prediction can be ‘uncertain’ in two different ways: (1)
the model’s predictions may have high variance, or (2) the model’s predictions
may be centered around a probability of interaction of $0.5$, where we are the
most unsure about whether this interaction exists. Improving the incorporation
of different forms of uncertainty in probabilistic interaction predictions seems
a necessary next step toward understanding what pairs of species we know the
least about, in order to prioritize sampling to provide the most new information
possible. 


# Positive associations in co-occurrence increase the false-negative-rate

The model above doesn't consider the possibility that there are positive or
negative associations which shift the probability of species cooccurrence away
from what is expected based on their relative abundances due to their
interaction [@Cazelles2016TheSpe]. However, here we demonstrate that the
probability of having a false-negative can be higher if there is some positive
association in the occurrence of species $A$ and $B$. If we denote the
probability that we observe the co-occurrence of two species $A$ and $B$ as
$P(AB)$ and if there is no association between the marginal probabilities of
observing $A$ and observing $B$, denoted $P(A)$ and $P(B)$ respectively, then
the probability of observing their co-occurrence is the product of the marginal
probabilities for each species, $P(AB) = P(A)P(B)$. In the other case where
there is some positive strength of association between observing both $A$ and
$B$ because this interaction is "important" for each species, then the
probability of observation both $A$ and $B$, $P(AB)$, is greater than $P(A)P(B)$
as $P(A)$ and $P(B)$ are not independent and instead are positively correlated,
i.e. $P(AB)> P(A)P(B)$. In this case, the probability of observing a single
false-negative in our naive model from @fig:geometric(A) is $p_{fn}= 1-P(AB)$,
which due to the above inequality implies $p_{fn}>1-P(A)P(B)$. This indicates an
increasingly greater probability of a false negative as the strength of
association gets stronger, $P(AB) \to P(AB) \gg P(A)P(B)$. However, this still
does not consider variation in species abundance in space and time
[@Poisot2015SpeWhy]. If positive or negative associations between species
structure variation in the distribution of $P(AB)$ across space/time, then the
spatial/temporal biases induced by data collection would further impact the
realized false-negative-rate, as the probability of false negative would not be
constant for each pair of species across sites.

To test for these positive associations in data we scoured Mangal for datasets
with many spatial or temporal replicates of the same system, which led the the
resulting seven datasets set in figure @fig:mangal. For each dataset,
we compute the marginal probability $P(A)$ of occurrence of each species $A$
across all networks in the dataset. For each pair of interacting species $A$ and
$B$, we then compute and compare the probability of co-occurrence if each
species occurs independently, $P(A)P(B)$, to the empirical joint probability of
co-occurrence, $P(AB)$. Following our analysis above, if $P(AB)$ is greater than
$P(A)P(B)$, then we expect our neutral estimates of the FNR above to
underestimate the realized FNR. In @fig:mangal, we see the difference between
$P(AB)$ and $P(A)P(B)$ for the seven suitable datasets with enough
spatio-temporal replicates and a shared taxonomic backbone (meaning all
individual networks use common species identifiers) found on Mangal to perform
this analysis. Further details about each dataset are reported in @tbl:id.

In each of these datasets, the joint probability of co-occurrence $P(AB)$ is
decisively greater than our expectation if species co-occur in proportion to
their relative abundance $P(A)P(B)$. This suggests that there may not be as many
“neutrally forbidden links” [@Canard2012EmeStr] as we might think, and that the
reason we do not have records of interactions between rare species is probably
due to observation error. This has serious ramifications for the widely observed
property of nestedness seen in bipartite networks
[@Bascompte2007PlaMut]---perhaps the reason we have lots of observations between
generalists is because they are more abundant, and this is particularly relevant
as we have strong evidence that generalism drives abundance [@Song2022GenDri],
not vice-versa.


![The difference between joint-probability of co-occurrence ($P(AB)$) and
expected probability of co-occurrence under independence ($P(A)P(B)$) for
interacting species for each dataset. The red-dashed line indicates 0 (no
association). Each histogram represents a density, meaning the area of the
entire curve sums to 1. The continuous density estimate (computed using local
smoothing) is shown in grey. The p-value on each plot is the result of a
one-sided t-test comparing the mean of each distribution to
0.](./figures/fig2.png){#fig:mangal}

|           Network          |      Type     |  $N$ |  $S$  |   $C$   |   $\bar{S}$   |   $\bar{C}$   |   $\bar{\beta}_{OS}$  |   $\bar{\beta}_{WN}$  |
|:--------------------------:|:-------------:|:---:|:---:|:-----:|:-----:|:-----:|:-----:|:-----:|
| @Kopelke2017FooStr   | Food Web      | 100 |  98 | 0.037 |  7.87 | 0.142 | 1.383 | 1.972 |
| @Thompson2000ResSol | Food Web      |  18 | 566 | 0.014 | 80.67 | 0.049 | 1.617 | 1.594 |
| @Havens1992ScaStr            | Food Web      |  50 | 188 | 0.065 | 33.58 | 0.099 | 1.468 | 1.881 |
| @Ponisio2017OppAtt      | Pollinator    | 100 | 226 | 0.079 |  23.0 | 0.056 | 1.436 | 1.870 |
| @Hadfield2014TalTwo     | Host-Parasite |  51 | 327 | 0.085 | 32.71 | 0.337 | 1.477 | 1.952 |
| @Closs1994SpaTem        | Food Web      |  12 |  61 |  0.14 | 29.09 | 0.080 | 1.736 | 1.864 |
| @CaraDonna2017IntRew    | Pollinator    |  86 | 122 |  0.18 | 21.42 | 0.312 | 1.527 | 1.907 |
Table: The datasets used in the above analysis (Fig 2). The table reports the type of each dataset, the total number of networks in each dataset $(N)$, the total species richness in each dataset $(S)$, the connectance of each metaweb (all interactions across the entire spatial-temporal extent) $(C)$, the mean species richness across each local network $\bar{S}$, the mean connectance of each local network $\bar{C}$, the mean $\beta$-diversity among overlapping species across all pairs of network species ($\bar{\beta}_{OS}$), and the mean $\beta$-diversity among all species in the metaweb ($\bar{\beta}_{WN}$). Both metrics are computed using KGL $\beta$-diversity [@Koleff2003MeaBet] {#tbl:id}


# The impact of false-negatives on network properties and prediction

Here, we assess the effect of false-negatives on our ability to make predictions
about interactions, as well as their effect on network structure. The prevalence
of false-negatives in data is the catalyst for interaction prediction in the
first place, and as a result methods have been proposed to counteract this bias
[@Stock2017LinFil; @Poisot2022NetEmb]. However, it is feasible that the FNR in a
given dataset is so high that it could induce too much noise for an interaction
prediction model to detect the signal of possible interaction between species.

To test this we use the dataset from @Hadfield2014TalTwo that describes
host-parasite interaction networks sampled across 51 sites, and the same method
as @Strydom2021RoaPre to extract latent features for each species in this
dataset based on applying PCA to the co-occurrence matrix. We then predict a
metaweb (equivalent to predicting true or false for an interaction between each
species pair, effectively a binary classification problem) from these
species-level features using four candidate models for binary
classification---three often used machine-learning (ML) methods (Boosted
Regression Tree (BRT), Random Forest (RF), Decision Tree (DT)), and one naive
model from classic statistics (Logistic Regression (LR)). Each of the ML models
are bootstrap aggregated (or bagged) with 100 replicates each. We partition the
data into 80-20 training-test split, and then seed the training data with false
negatives at varying rates, but crucially do nothing to the test data. We fit
all of these models using MLJ.jl, a high-level Julia framework for a
wide-variety of ML models [@Blaom2020MljJul]. We evaluate the efficacy of these
models using two common measures of binary classifier performance: the area
under the receiver-operator curve (ROC-AUC) and the area under the
precision-recall curve (PR-AUC), for more details see @Poisot2022GuiPre. Here,
PR-AUC is slightly more relevant as it is a better indicator of prediction of
false-negatives. The results of these simulations are shown in
@fig:addedfnr(A&B).

![**(A)** The area-under the receiver-operator curve (ROC-AUC) and **(B)** The area-under the precision-recall curve (PR-AUC; right) for each different predictive model (colors/shapes) across a spectrum of the proportion of added false-negatives (x-axis). **(C)** The mean trophic-level of all species in a network generated with the niche model across different species richnesses (colors). For each value of the FNR, the mean trophic level was computed across 50 replicates. The shaded region for each line is one standard-deviation across those replicates.](./figures/fig3.png){#fig:addedfnr}

One interesting result seen in @fig:addedfnr(A&B) is that the ROC-AUC value does
not approach random in the same way the PR-AUC curve does as we increase the
added FNR. The reason for this is that ROC-AUC is fundamentally not as useful a
metric in assessing predictive capacity as PR-AUC. As we keep adding more
false-negatives, the network eventually becomes a zeros matrix, and these models
can still learn to predict “no-interaction” for all possible species pairs,
which does far better than random guessing (ROC-AUC = 0.5) in terms of the false
positive rate (one of the components of ROC-AUC). This highlights a more broad
issue of label class imbalance, meaning there are far more non-interactions than
interactions in data. A full treatment of the importance of class-balance is
outside the scope of this paper, but is explored in-depth in @Poisot2022GuiPre.
Further we see, if anything, gradual decline in the performance of the model
until we reach very high FNR levels (i.e. $p_{fn} > 0.7$). This is consistent
with other recent work [@Gupta2023HowMan], although it must be considered that
the empircal data on which these models are trained already are almost certain
to already contain false-negatives. 

Although these ML models are surprisingly performant at link prediction given
their simplicity, there have been several major developments in applying
deep-learning methods to many tasks in network inference and prediction---namely
graph-representation learning (GRL, @Khoshraftar2022SurGra) and graph
convolutional networks [@Zhang2019GraCon]. At this time, these advances can not
yet be applied to ecological networks because they require far more data than we
currently have. We already have lots of features that could be used as inputs
into these models (i.e. species level data about occurrence, genomes, abundance,
etc.), but our network datasets barely get into the hundreds of local networks
sampled across space and time (@tbl:id). Once we start to get into the
thousands, these models will become more useful, but this can only be done with
systematic monitoring of interactions. This again highlights the need to
optimize our sampling effort to maximize the amount of information contained in
our data given the expense of sampling interactions.

We also consider how the FNR affects network properties. In @fig:addedfnr(C) we
see the mean trophic level across networks simulated using the niche model (as
above), across a spectrum of FNR values. In addition to the clear dependence on
richness, we see that mean trophic level, despite varying widely between niche
model simulations, tends to be relatively robust to false-negatives and does not
deviate widely from the true value until very large FNRs.
This is not entirely unsurprising. Removing links randomly from a food-web is
effectively the inverse problem of the emergence of a giant component (more than
half of the nodes are in a connected network) in random graphs (see
@Li2021PerCom for a thorough review). The primary difference being that we are
removing edges, not adding them, and thus we are witnessing the dissolution of a
giant component, rather than the emergence of one. Further applications of
percolation theory [@Li2021PerCom] to the topology of sampled ecological
networks could improve our understanding of how false-negatives impact the
inferences about the structure and dynamics on these networks.

# Discussion

Species interactions enable the persistence and functioning of ecosystems, but
our understanding of interactions is limited due to the intrinsic difficulty of
sampling them. Here we have provided a null model for the expected number of
false-negatives in an interaction dataset. We demonstrated that we expect many
false-negatives in species interaction datasets purely due to the intrinsic
variation of abundances within a community. We also, for the first time to our
knowledge, measured the strength of association between co-occurrence and
interactions [@Cazelles2016TheSpe] across many empirical systems, and found that
these positive associations are both very common, and showed algebraically that
they increase the realized FNR. We have also shown that false-negatives could
further impact our ability to both predict interactions and infer properties of
the networks, which highlights the need for further research into methods for
correcting this bias in existing data.

A better understanding of how false-negatives impact species interaction data is
a practical necessity---both for inference of network structure and dynamics,
but also for prediction of interactions by using species level information.
False-negatives could pose a problem for many forms of inference in network
ecology. For example, inferring the  dynamic stability of a network could be
prone to error if the observed network is not sampled "enough". What exactly
"enough" means is then specific to the application, and should be assessed via
methods like those here when designing samples. Further, predictions about
network rewiring [@Thompson2017DisGov] due to range shifts in response to
climate change could be error-prone without accounting for interactions that
have not been observed but that still may become climatically infeasible. As is
evident from @fig:geometric(A), we can never guarantee there are no
false-negatives in data. In recent years, there has been interest toward
explicitly accounting for false-negatives in models [@Stock2017LinFil;
@Young2021RecPla], and a predictive approach to networks---rather than expecting
our samples to fully capture all interactions [@Strydom2021RoaPre]. As a result,
better models for predicting interactions are needed for interaction networks.
This includes explicitly accounting for observation error
[@Johnson2021BayEst]---certain classes of models have been used to reflect
hidden states which account for detection error in occupancy modeling
[@Joseph2020NeuHie], and could be integrated in the predictive models of
interactions in the future.

This work has several practical consequences for the design of surveys for
species' interactions. Simulating the process of observation could be a powerful tool for
estimating the sampling effort required by a study that takes relative abundance
into account, and provides a null baseline for expected FNR. It is necessary to
take the size of the species pool into account when deciding how many total
samples is sufficient for an “acceptable” FNR (@fig:geometric(C & D)). Further
the spatial and temporal turnover of interactions means any approach to sampling
prioritization must be spatiotemporal. We demonstrated earlier that observed
negatives outside of the range of both species aren’t informative, and therefore
using species distribution models could aid in this spatial prioritization of
sampling sites.

We also should address the impact of false-negatives on the inference of process
and causality in community ecology. We demonstrated that in model food webs,
false-negatives do not impact the measure of total trophic levels until very
high FNR (figure @fig:addedfnr(C)), although we cannot generalize this further
to other properties. This has immediate practical concern for how we design what
taxa to sample---does it matter if the sampled network is fully connected? It
has been shown that the stability of subnetworks can be used to infer the
stability of the metaweb paper beyond a threshold of samples [@Song2022RapMon].
But does this extend to other network properties? And how can we be sure we are
at the threshold at which we can be confident our sample characterizes the whole
system? We suggest that modeling observation error as we have done here can
address these questions and aid in the design of samples of species
interactions. To try to survey to avoid all false-negatives is a fool's errand.
Species ranges overlap to form mosaics, which themselves are often changing in
time. Communities and networks don't end in space, and the interactions that
connect species on the 'periphery' of a given network to species outside the
spatial extent of a given sample will inevitably appear as false-negatives in
practical samples. The goal should instead be to sample a system enough to have
a statistically robust estimate of the current state and empirical change over
time of an ecological community at a given spatial extent and temporal
resolution, and to determine what the sampling effort required should be prior to
sampling.

# Acknowledgements

AG & MDC acknowledge the support of the Liber Ero Chair for Biodiversity
conservation and NSERC.  

# References
