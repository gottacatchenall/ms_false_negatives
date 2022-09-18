---
bibliography: [references.bib]
---

# Introduction

Community ecology is fundamentaly about understanding species interactions.
Understanding how species interact underlies many  question in
evolutionary biology, community ecology, and is an increasing imperative to mitigate the anthropogoenic
change on  Earth's biodiversity [@Makiola2020KeyQue; @Jordano2016ChaEco] and to predict
potential spillover of zoonotic disease [@Becker2021OptPre]. Over the past
couple decades biodiversity data has become increasingly
available---remote-sensing has enabled collection of data on spatial scales and
resolutions previously unimaginable, and improved in-situ sensing
[@Stephenson2020TecAdv] and adoption of open data practices [@Kenall2014OpeFut]
have substantially amount of data available to ecologists. Still, widespread
data about species _interactions_ remains elusive  [@TIM_DATASET_PAPER]. Often,
observing an interaction between two species requires human sampling, typically
because remote sampling methods can only detect co-occurrence, and this itself
is not necessarily indicative of interaction [@Blanchet2020CooNot]. This
constraint induces biases on species interaction data subject to the spatial and
temporal scales that humans can feasibly sample.

_Sampling effort_ and its impact on the resulting data collected from ecosystems
has encouraged a long history of discourse. The recorded number of species in a
sample depends on the total number of observations [@Willott2001SpeAcc;
@Walther1995SamEff], as do estimates of population abundance
[@Griffiths1998SamEff]. This has motivated more quantitatively robust approaches
to account for error in sampling data in many contexts: to determine if a given
species is extinct [@Boakes2015InfSpe], to determine sampling design
[@Moore2016OptEco], and to measure species richness across large scales [@Carlson2020WhaWou].
In the context of interactions, the initial concern was the compounding effects
of limited sampling effort combined with the amalgamation of data (across both
study sites and across taxonomic scales) could lead any empirical set of
observations to inadequately reflect the reality of how species interact
[@Paine1988RoaMap]. @Martinez1999EffSam showed that in a plant-endophyte trophic
network, network connectance is robust to sampling effort, but this done in the
context of a system for which observation of 62,000 total interactions derived
from 164,000 plant-stems was feasible. In some systems (e.g. megafauna
food-webs) this many observations is either impractical or infeasible due to the
absolute abundance of the species in question.

Because we cannot feasibly observe all (or even many) of the interactions that occur in
nature, our samples end up capturing only a small fraction of those
interactions. This means we can be reasonably confident two species actually
interact if we have a record of it, but not at all confident that two species
_do not_ interact if we have no record of those species observed together. In
other words, it is difficult to distinguish true-negatives (two species _never_ interact)
from _false-negatives_ (two species interact in some capacity, but we have not
observed it). This is then amplified as the interaction data we have is
geographically toward the usual suspects [@Poisot2021GloKno], This noise in data
has practical consequences for answering questions about species interactions
[@deAguiar2019RevBia]---these false-negatives could go on to effect the
inferences we make about network properties and relations among species, and our
predictions about how species will interact in the future.

This is compounded by semantic confusion about the definition of "interaction".
Here distinguish between: a species _occurring_, a species being _observed
occurring_, two species being observed _co-occurring_, and two species being
observed _interacting_ (@fig:taxonomy). In this manuscript, we refer to species
either as "interacting"---two species co-occur (and, at least sometimes,
interact)---or "not-interacting" (two species that, regardless of whether they
co-occur, neither exhibits any meaningful effect on the biomass of the other).
In @fig:taxonomy we see that, under our definition, observing two species
co-occurring is a prerequisite for observing an interaction between two species.


But species are not observed with equal probability but instead in proportion to
their relative biomass---you are much more likely to observe a species of high
relative abundance than one of very low relative abundance [@Poisot2015SpeWhy].
This assumes that there are no associations in species co-occurrence due to an
interaction (perhaps because this interaction is "important" for both species)
[@Cazelles2016TheSpe], but here we show increasing strength of associations
leads to increasing probability of false-negatives in interaction data. Further
observed co-occurrence is often equated with meaningful interaction strength,
but this is not necessarily the case [@Blanchet2020CooNot; @Strydom2021RoaPre].
Bears and salmon _interact_---a bear and the microbes in the soil of a dens
interact, but less so.

Here, we show that the probability of observing a actual "non-interaction" between
species depends on sampling effort, and suggest that surveys of species
interactions can benefit from simulation modeling of detection probability
[@Jordano2016SamNet]. We demonstrate that the realized false-negative rate of
interactions is directly related the relative abundance of a particular species,
relationship between total sampling effort (the total count of all individuals
of all species seen) and false-negative rate. questions we pose and attempt to
answer are: 1) How many times do you have to observe a non-interaction between
two species to be confident in saying that is a true negative? 2) How "wrong"
are the measurements of network structure as a function of false-negative
probability? and lastly 3) How do false-negatives impact our ability to make
reliable predictions about interactions? We show that positive associations in
co-occurrence data can increase realized probability of false negatives, and
demonstrate these positive associations are present in two spatially-replicated
systems. We conclude by suggesting that simulation of sampling effort and
species occurrence can and should be used to help design surveys of species
diversity [@Moore2016OptEco], and by advocating use of null models like those
presented here as a tool for guiding design of surveys of species interactions,
and for modeling detection error in predictive ecological models.

# How many observations of a non-interaction do we need to classify it as a true negative?

To answer the titular question of this section, we present a naive model of
interaction detection: we assume that every interacting pair of species is
incorrectly observed as a not-interacting with an independent and fixed
probability, which we denote $p_{fn}$ and subsequently refer to as the
False-Negative Rate (FNR). If we observe the same species not-interacting $N$
times, then the probability of a true-negative (denoted $p_{tn}$) is given by
$p_{tn} = 1 - (p_{fn})^N$. This relation (callend the geometric distribution, a
special case of the negative-binomial distribution) is shown in
@fig:negativebinom for varying values of the false negative rate $p_{fn}$. This
illustrates a fundamental link between our ability to reliably say an
interaction doesn't exist---$p_{tn}$---and the number of times we have observed
a given species. In addition, note that there also is no non-zero $p_{fn}$ for
which we can ever _prove_ that an interaction does not exist---no matter how
many observations of non-interaction $N$ we have, $p_{tn} < 1$.

From @fig:fig1 (A) (and general intuition) it is clear that the more times
we see two species _occurring_, but _not_ interacting, the more likely the
interaction is a true negative. But how does one decide what this threshold of
number of observations should be when planning to sample a given system? If
false-negative rates presented in @fig:fig1 seem unrealistically high,
consider that species are not observed independent of their relative abundance.
In the next section we demonstrate that distribution of abundance in ecosystems
can lead to realized values of $p_{fn}$ similar to those in @fig:negativebinom
for species with low relative abundance, simply as a artifact of sampling
effort.

## False-negatives as a product of relative abundance

Here we show the realized false-negative rate of species interactions changes
drastically with sampling effort, largely due to the intrinsic variation of
abundances within a community. We do this by simulating the process of
observation of species interactions, applied both to 243 empirical food webs
from the Mangal database [@Banville2021ManJl] as well as random food-webs
generated using the niche model [@Williams2000SimRul]. Our neutral model of
observation assumes each observed species is drawn from the distribution of
those species' abundances at that place and time. Although there is no shortage
of debate as to the processes the govern this distribution of abundances within
a community, this abundance distribution can be reasonably-well described by a
log-normal distribution [@Volkov2003NeuThe]  (Note that in addition to the
log-normal distribution, we also tested the case where the abundance
distribution is derived from power-law scaling $Z^{(T_i-1)}$ where $T_i$ is the
trophic level of species $i$  and $Z$ is a scaling coefficient.
[@Savage2004EffBod], which yields the same qualitative behavior, _supplement figure
1_). The practical consequence of abundance distributions spanning many order of
magnitude is seeing two "rare" species interacting requires two low
probability events: observing two rare species _at the same
time_.

To simulate the process of observation, for an ecological network $A$ with $S$
species, we sample abundances for each species from a standard-log-normal
distribution. For each true interaction in $A$ (i.e. $A_{ij} = 1$) we estimate
the probability of observing both species $i$ and $j$ at given place and time by
simulating $n$ observations of individuals, where the species of the individual
observed at the $1,2,\dots,n$-th observation is drawn from the generated
log-normal distribution of abundances. For each pair of species $(i,j)$, if both
$i$ and $j$ are observed within the $n$ observations, the interaction is tallied
as a true positive if $A_{ij}=1$ and a false positive otherwise. Similarly, if
only one of $i$ and $j$ are observed---_but not both_---in these $n$
observations, but $A_{ij}=1$, this is counted as a false-negative, and a
true-negative otherwise.

In @fig:fig1 (C) and (D) we see this model of observation applied to networks
generated using the niche model [@Williams2000SimRul] across varying levels of
species richness, and in (b) applied to 243 food-webs from the Mangal database.
For all niche model simulations in this manuscript, for a given number of
species $S$ the number of interactions is drawn from the flexible-links model
fit to Mangal data [@MacDonald2020RevLin], effectively drawing the number of
interactions $L$ for a random niche model food-web as $L \sim
\text{BetaBinomial}(S^2-S+1, \mu \phi, (1-\mu)\phi)$, where the MAP estimate of
($\mu$, $\phi$) applied to Mangal data from @MacDonald2020RevLin is $(\mu =
0.086, \phi =24.3)$. All simulations were done with 500 independent replicates
per unique number of observations $n$. All analyses presented here are done in
Julia v1.6 [@Bezanson2015JulFre] using both EcologicalNetworks.jl v0.5 and
Mangal.jl v0.4 [@Banville2021ManJl; ZENODO link TODO]. Note that the empirical
data also is, due to the phenomena described here, very likely to _already_ have
many false negatives, which is why we are interested in prediction of networks
in the first place---we'll revisit this in the final section.

![A) The probability an observed interaction is a true negative (y-axis) given how many times it has been sampled as a non-interaction (x-axis). Each color reflects a different value of $p_{fn}$, the false-negative rate (FNR). this is effectively the cdf of the negative-binomial distribution with $r=1$. (B) The expected needed observations of all individuals of all species (y-axis) required to obtain a goal number of observations (colors) of a particular species, and a function of the relative abundance of that focal species (x-axis). (C) and (D): False negative rate (y-axis) as a function of total sampling effort (x-axis) and network size, computed using the method described above. For 500 independent draws from the niche model Williams2000SimRul at varying levels of species richness (colors) with connectance drawn according to the flexible-links model MacDonald2020RevLin as described in the main text. For each draw from the niche model, 200 sets of 1500 observations are simulated, for which each the mean false negative rate at each observation-step is computed. Means denoted with points, with $1\sigma$ in the first shade and $2\sigma$ in the second. B: empirical food webs from Mangal database in teal, applied to the same process as the A. The outlier on panel B is a 714 species food-web](./figures/fig1.png){#fig:fig1}

In panel (c) of @fig:totalobs, we show the expected number of total observations
needed to obtain a "goal" number of observations (colors) of a particular
"focal" species. As an example, if we hypothesize that $A$ and $B$ do not
interact, and we want to see species $A$ and $B$ both co-occurring and
not-interacting 10 times to be confident this is a negative (a la
@fig:negativebinom), then we need an expected 10,000 observations of all species
if the relative abundance of $A$ is 0.00125.

Empirical data on interactions are subject to the practical limitations of
funding and human-work hours, and therefore existing data tend to fall on the
order on 100s or 1000s observations of individuals per site [@Resasco2021PlaPol;
@Schwarz2020TemSca; @Nielsen2007EcoNet]. Clear aggregation of this data has
proven difficult to find and a meta-analysis of network data and sampling effort
seems both pertinent and necessary, in addition to the effects of aggregation of
interactions across taxonomic scales [@Giacomuzzo2021FooWeb;
@Gauzens2013FooAgg]. Further, from @fig:totalobs it is evident that the number
of species considered in a study is inseparable from the false-negative rate in
that study, and this effect should be taken into account when designing samples
of ecological networks in the future.

We conclude this section by advocating for the use of neutral models similar to
above to generate expectations about the number of false-negatives in a data set
of a given size. This could prove fruitful both for designing surveys of
interactions [@Canard2012EmeStr],  but also because we may want to incorporate
models of  observation error into predictive models [@Joseph2020NeuHie].
Additionally, one must consider the context for sampling---is the goal to detect
a particular species $A$ (as in @fig:totalobs (c)), or to get a representative
sample of interactions across the species pool? This argument is well-considered
when sampling species [@Willott2001SpeAcc], but has not yet been internalized
for designing samples of communities.

## Positive associations can increase the probability of false-negatives

This model above doesn't consider the possibility that there are positive or
negative associations which shift the probability of observing two species
together due to their interaction [@Cazelles2016TheSpe]. However, here we
demonstrate that the probability of observing a false negative can be _higher_
if there is some positive association between occurrence of species $A$ and $B$.

If we denote the probability that we observe an interaction we know exists
between $A$ and $B$ as $P(AB)$, and if there is _no_ association between the
marginal probabilities of observing $A$ and observing $B$, denoted $P(A)$ and
$P(B)$ respectively, then the probability of observing the interaction $P(AB) =
P(A)P(B)$. In the other case where there _is_ some positive strength of
association between observing both $A$ and $B$ because this interaction is
"important" for each species, then the probability of observation both $A$ and
$B$, $P(AB)$, is greater than $P(A)P(B)$ as $P(A)$ and $P(B)$ are not
independent and instead are positively correlated, _i.e._ $P(AB) > P(A)P(B)$. In
this case, the probability of observing a false negative in our naive model from
@fig:negativebinom is $p_{fn} = 1 - P(AB)$ which due to the above inequality
implies $p_{fn} \geq 1 - P(A)P(B)$ which indicates increasingly greater
probability of a false negative as $P(AB) \to P(AB) \gg P(A)P(B)$.

However this does not consider variation in species abundance in space and time,
[@Poisot2015SpeWhy]. If positive or negative associations between species
structure variation in the distribution of $P(AB)$ across space/time, then the
spatial/temporal biases induced by data collection would further impact the
realized false negative rate, as the probability of false negative would not be
constant for each pair of species across sites. To test for this association
empirical data, we use two datasets: a set of host-parasite interactions sampled
across 51 sites with 327 total taxa  [@Hadfield2014TalTwo] and a set of 18 New
Zealand freshwater stream food webs with 566 total taxa [@Thompson2000ResSol].
We simply compute the empirical marginal distribution of species occurrence, and
compare the product of the marginals, $P(A)P(B)$, to the empirical joint
distribution $P(AB)$.

![Top: Hadfield, Bottom: NZ Stream Foodwebs. Effectively a version of
@Cazelles2016TheSpe figure 1 panel
A. Both distributions have $\mu \neq 0$ with $p < 10^{-50}$](./figures/positiveassociations.png){#fig:associations}

In @fig:associations, both host-parasite system (top) and food-web (bottom)
exhibit these positive associations. There is no reason to expect the strength
of this association to be the same in different systems. At the moment,
computing this metric for all of the networks in the Mangal database proves
challenging as most data sets use different taxonmic identifiers, often at
different resolutions. These particular datasets [@Hadfield2014TalTwo;
@Thompson2000ResSol] were usable because they already have been sorted to have a
fixed taxonomic backbone (as part of EcologicalNetworks.jl
[@Banville2021ManJl]). Applying this in bulk to Mangal food-webs presents the
difficulty of resolving different taxon identifiers across spatial samples of
species with to different resolutions, which is why we can't simply apply this
to the whole Mangal database---this highlights a general problem of resolving
taxonomic indentifiers which use different names and different resolutions in
different ecological datasets, which is a problem that needs to be addressed for
computational approaches to scale up to the world of big-ecological-data we hope
to build, although this is a task that may be aided via
natural-language-processing methods.


# The impact of false-negatives on network analysis and prediction

We now transition toward assessing the effects of false negatives in data on
the properties of the networks which we derive from this interaction data, and
their effect on models for predicting interactions in the future.

## Effects of false-negatives on network properties

Here we simulate the process of observation with error to generate synthetic
data with a known proportion of false negatives, and compare the computed
network properties of the original "true" network to the computed properties of
the "observed" network with added false-negatives. In @fig:properties we see the
mean-squared error of connectance, mean degree-centrality, and spectral radius,
computed across 2000, 2000, and 300 replicates respectively at each value of the
false negative rate $p_{fn}$. All replicates use random food-webs simulated
using the niche model [@Williams2000SimRul] with $100$ species and connectance
drawn from the flexible-links model [@MacDonald2020RevLin] as before.

![The mean-squared error (y-axis) of various network properties (different
colors) across various simulated false-negative rates (x-axis). Means denoted
with points, with $1\sigma$ in the first shade and $2\sigma$ in the
second.](./figures/props_specrad.png){#fig:properties}

We consider three properties: connectance, mean-degree-centrality, and spectral
radius, indicative of local, meso, and global structure. Connectance is
effectively a node-level property, a proxy for the degree distribution.
Degree-centrality captures a different aspect of network structure than
connectance, more indicative of meso-level properties that describe local
'regions' of nodes interact. Spectral radius (equivalent to the magnitude of the
largest eigenvalue of $A$) is a measure of global structure, and demonstrates
the most variability in response to false-negatives. For example, if a
false-negative splits a metaweb into two components, spectral-radius becomes the
largest eigenvalue of each of those two components. Also note that the form of
this error function varies little as species richness changes (_supplemental
figure 2_). Practically, @fig:properties shows us that different scales of
measuring network structure vary in their response to false
negatives---connectance responds roughly linearly to false negatives, whereas
mean-degree-centrality decisively does not. This implies that false-negatives
adversely could effect indirect interactions [@Williams2002TwoDeg].

## Effects of false negatives on ability to make predictions

Here, we assess the effect of false negatives in data on our ability to make
predictions about interactions. The prevalence of false-negatives in data is the
catalyst for interaction prediction in the first place, and as a result methods
have been proposed to counteract this bias [@Poisot2021ImpMam;
@Stock2017LinFil]. However, it is feasible this could induce too much noise for
a interaction prediction model to detect the signal of interaction chance from
to the latent properties of each species derived from the empirical network if
the number of false-negatives in a dataset becomes too overwhelming.

To test this, we use the same predictive model and dataset as in
@Strydom2021RoaPre to predict a metaweb from various empirical slices of the
species pool observed across space. This dataset from @Hadfield2014TalTwo
describes host-parasite interaction networks sampled across 51 sites. We
partition the data into 80-20 training-test split, and then seed the training
data with false negatives varying rates, but crucially do nothing to the test
data. We use the same model, a neural-network with 3 feed-forward layers to
predict outputs based on features extracted from co-occurence (see
@Strydom2021RoaPre for more details). The single modification we make to the
model is not enforcing a number of positives in the training data as this
constraint is  eventually impossible for increasing FNR. In @fig:rocpr, we show
receiving-operating-characteristic (ROC) and precision-recall (PR) curves for
the model with varying levels of synthetic false-negatives added to the data.

![Receiver-operating-characteristic (left) and precision-recall (right) curves
for the model on varying levels of false-negatives in the data (colors). For
each value of FNR, we run 30 random training/test splits on 80/20 percent of the
data. Replica of figure 1 in
@Strydom2021RoaPre](./figures/rocpr_falsenegatives.png){#fig:rocpr}

Interestingly, the performance of the model from @Strydom2021RoaPre changes
little with many added false-negatives, which is good evidence in favor
neural-networks as a class of model for interaction detection. Again, similar to
our caveat in the previous section, this dazta is _already_ likely to have many
false-negatives, so the effects of adding more as we do in this illustration
might be mitigated because there are already non-simulated false-negatives in
the original data which impact the models performance, even in the $p_{fn} = 0$
case.







We conclude be proposing that simulating the effects of false negatives in this
way can serve as an additional validation tool when aiming to detect structural
properties of networks using generative null models [@Connor2017UsiNul], or when
evaluating the robustness of a predictive model.   

# Discussion

Here, we have demonstrated that we expect false-negatives in species interaction
datasets purely due to the distribution of abundances within a community.
Positive associations between species occurrence [@Cazelles2016TheSpe] can
increase the realized false-negative rate if the sampling effort is limited, and
we have presented evidence of this non-random structure of co-occurrence in two
sets of spatially-replicated ecological network samples. We have also shown that
false-negatives can cause varying responses in our measurements of network
properties and further could impact our ability to reliably predict
interactions, which highlights the need for further research into methods for
correcting this bias in existing data [@Stock2017LinFil]. A brief caveat here is
that we do not consider the rate of false-positives---in large part
false-positives can be explained by misidentification of species, although this
could be a relevant consideration in some cases.

What does the future hold for this research? A better understanding of how
false-negatives impact our analyses and prediction of ecological networks is a
practical necessity. False-negatives could pose a problem for many forms of
inference in network ecology. For example, if we aim to measure structural or
dynamic stability of a network, or to infer indirect interactions
[@Williams2002TwoDeg], these estimates could be prone to error if the observed
network is not sampled "enough". What exactly "enough" means is then specific to
the application, and should be assessed via methods like those here when
designing samples. Further, predictions about network rewiring
[@Thompson2017DisGov] due to a changing climate could be error-prone without
accounting for interactions that have not been observed but that still may
become climatically infeasible.  

This highlights the need for a quantitatively robust approach samples design:
for interactions [@Jordano2016SamNet] and otherwise [@Carlson2020WhaWou]. The
primary takeaway is that when planning the sampling effort across sites, it is
necessary to take both the size of the species pool into account. Further,
simulating the process of observation could be a powerful tool for planing study
design which takes relative abundance into account, and provide a null baseline
for detection of interaction strength. A model similar to that here can and
should be used to provide a neutral expectation of true-negative probability
given a number of observations of individuals at a given place and time.

As we derive from @fig:negativebinom, we can never guarantee there are no
false-negatives in data. In recent years, there has been interest toward
explicitly accounting for false-negatives in models [@Young2021RecPla;
@Stock2017LinFil], and toward a predictive approach toward interactions
---rather than expect that our samples can fully capture all interactions, we
know that some interactions between species will not be observed due to finite
sampling capacity, and instead we must impute  the true metaweb of interactions
given a set of samples [@Strydom2021RoaPre]. As a result, better predictive
approaches are needed for interaction networks [@Strydom2021RoaPre], and
building models that explicitly account for observation error is a necessary
step forward for predictive ecological models [@Young2021RecPla;
@Johnson2021BayEst]. Neural networks, like the one used to predict interactions
in the above section, have been used to reflect hidden states which account for
detection error in occupancy modeling [@Joseph2020NeuHie], and could be
integrated in the predictive models of the future.

A better conceptual framework for designing surveys and monitoring networks, and
incorporating sequential observations over time is clearly needed
[@Carlson2020WhaWou], combined with a meta-analysis of sampling effort and
taxonomic resolution in existing data. Incorporating a better understanding of
sampling effects and bias on both the future design of biodiversity monitoring
systems, and the predictive models we wish to apply to this data, is imperative
in making actionable predictions about the future of ecological interactions on
our planet.


# References
