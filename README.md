---
bibliography: [references.bib]
---

# Introduction

Species interactions drive many processes in evolutionary biology and community
ecology. A better understanding of interactions among species is an imperative
to both mitigate the potentially harmful impacts of anthropogenic change on
Earth's biodiversity [@Makiola2020KeyQue] and to predict zoonotic spillover of
disease to prevent future pandemics [@Becker2021OptPre]. However, meeting these
challenges is difficult because interactions are intrinsically hard to sample
[@Jordano2016SamNet]. Over the past few decades biodiversity data has become
increasingly available---for example, remote-sensing has enabled collection of
data on spatial scales and resolutions previously unimaginable
[@Stephenson2020TecAdv], while the adoption of open data practices
[@Kenall2014OpeFut] have substantially increased the amount of data available to
ecologists. Still, widespread data about species interactions remains elusive
[@Poisot2021GloKno]. Observing an interaction between two species often requires
human observation because remote sampling methods can primarily detect
co-occurrence (but see @Niedballa2019AssAna), and this itself is not necessarily
indicative of an interaction [@Blanchet2020CooNot]. This constraint induces
biases on species interaction data subject to the spatial and temporal scales
that current observation methods can feasibly sample. This is further compounded
by semantic confusion around the word “interaction”---for example one might
consider competition a type of species interaction, even though it is marked by
a lack of co-occurrence between species, unlike other types of interactions,
like trophism or pollination, which require both species to be together at the
same place and time. We define interaction in the latter sense, where two
species have fitness consequences on one-another if they are in the sample place
at the same time. In addition, here we only consider direct (not higher-order)
interactions.

The importance of _sampling effort_ and its impact on resulting ecological data
has produced a rich body of literature. The recorded number of species in a
dataset or sample depends on the total number of observations
[@Willott2001SpeAcc; @Walther1995SamEff]---as do estimates of population
abundance [@Griffiths1998SamEff]---in addition to spatial coverage and species
detectability. This has motivated more quantitatively robust approaches to
account for error in sampling data in many contexts: to determine if a given
species is extinct [@Boakes2015InfSpe], to determine sampling design
[@Moore2016OptEco],and to measure species richness across large scales
[@Carlson2020WhaWou]. In the context of interactions, an initial concern was the
compounding effects of limited sampling effort combined with the amalgamation of
data (across both study sites, time of year, and taxonomic scales) could lead
any empirical set of observations to inadequately reflect the reality of how
species interact [@Paine1988RoaMap] or the structure of the network as a whole
[@McLeod2021SamAsy]. @Martinez1999EffSam showed that in a plant-endophyte
trophic network, network connectance is robust to sampling effort, but this was
done in the context of a system for which observation of 62,000 total
interactions derived from 164,000 plant-stems was feasible. In some systems
(e.g. megafauna food-webs) this many observations is either impractical or
infeasible due to the absolute abundance of the species in question.

We cannot feasibly observe all (or even most) of the direct interactions that
occur in an ecosystem. This means we can be confident two species actually
interact if we have a record of it (given an estimate of species
misidentification probability), but not at all confident that a pair of species
_do not_ interact if we have no record of those species observed together. In
other words, it is difficult to distinguish _true-negatives_ (two species never
interact) from _false-negatives_ (two species interact sometimes, but we do not
have a record of it). For a concrete example of a false-negative in a food web,
see @fig:concept. Because even the most highly sampled systems will still
contain missing interactions, there is increasing interest in combining
species-level data (e.g. traits, abundance, range, phylogenetic relatedness,
etc.) to build models to predict interactions between species we haven't
observed together before [@Strydom2021RoaPre]. However, the noise of
false-negatives could impact the efficacy of our predictive models and have
practical consequences for answering questions about interactions
[@deAguiar2019RevBia]. This data constraint is amplified as the interaction data
we have is geographically biased toward the usual suspects [@Poisot2021GloKno].
We therefore need a systematic approach to assessing these biases in the
observation process and the consequences this has for our understanding of
interaction networks.

The intrinsic properties of ecological communities create several challenges for
sampling: first, species are not observed with equal probability---we are much
more likely to observe a species of high abundance than one of very low
abundance [@Poisot2015SpeWhy]. @Canard2012EmeStr presents a null model of
food-web structure where species encounter one-another directly in proportion to
each species’ relative-abundance. This assumes that there are no associations in
species co-occurrence due to an interaction (perhaps because this interaction is
"important" for both species; @Cazelles2016TheSpe), but in this paper we later
show increasing strength of associations leads to increasing probability of
false-negatives in interaction data, and that these positive associations are
rampant in existing network data. Second, observed co-occurrence is often
equated with meaningful interaction strength, but this is not necessarily the
case [@Blanchet2020CooNot]---a true “non-interaction" would require that neither
of two species, regardless of whether they co-occur, exhibit any meaningful
effect on the fitness of the other. So, although co-occurrence is not directly
indicative of an interaction, it is a precondition for an interaction. Therefore
observations of “non-interactions” between pairs of species that are outside of
the union of both species ranges do not provide any information about that
interaction, i.e. they should be excluded from consideration.

Here, we illustrate how our confidence that a pair of species we believe to not
interact highly depends on sampling effort, and suggest that surveys of species
interactions can benefit from simulation modeling of detection probability. We
demonstrate that the realized false-negative rate of interactions is directly
related to the relative abundance of the species involved, and demonstrate how
simulation can be used to produce a null estimate of the false-negative rate as
a function of total sampling effort (the total count of all individuals of all
species seen). We show that positive associations in co-occurrence data can
increase realized probability of false-negatives, and demonstrate these positive
associations are ubiquitous in network datasets. We conclude by recommending
that the simulation of sampling effort and species occurrence can and should be
used to help design surveys of species interaction diversity [@Moore2016OptEco],
and by advocating use of null models like those presented here as a tool for
both guiding design of surveys of species interactions and for modeling
detection error in predictive models.

![This conceptual example considers a sample of the trophic community of bears, wolves, salmon (pink fish), pike (yellow fish), berry trees, and aspen trees. The true metaweb (all realized interactions across the entire spatial extent) is shown on the left. In the center is what a hypothetical ecologist samples at each site. Notice that although bears are observed co-occurring with both salmon and pike, there was never a direct observation of bears eating pike, even though they actually do. Therefore, this interaction between bears and pike is a false negative.](./figures/concept.png){#fig:concept}

# How many observations of a non-interaction do we need to be confident it's a true negative?

We start with a naive model of interaction detection: we assume that every
interacting pair of species is incorrectly observed as not-interacting with an
independent and fixed probability, which we denote pfn and subsequently refer to
as the False-Negative Rate (FNR). If we observe the same species not-interacting
$N$ times, then the probability of a true-negative (denoted $p_{tn}$) is given by
$p_{tn}=1-(p_{fn})^N$. This relation (the probability-mass-function of geometric
distribution, a special case of the negative-binomial distribution) is shown in
Figure 1(A) for varying values of $p_{fn}$ and illustrates a fundamental link
between our ability to reliably say an interaction doesn't
exist---$p_{tn}$---and the number of times we have observed a given species. In
addition, note that also there is no non-zero $p_{fn}$ for which we can ever
prove that an interaction does not exist---no matter how many observations of
non-interactions $N$ we have, $p_{tn}<1$.

From @fig:geometric(A) it is clear that the more often we see two species
co-occurring, but not interacting, the more likely the interaction is a true
negative. This has several practical consequences: first it means negatives
taken outside the overlap of the range of each species aren’t informative.
Second, we can use this relation to compute the expected number of total
observations needed to obtain a "goal" number of observations of a particular
pair of species (@fig:geometric(B)). As an example, if we hypothesize that $A$
and $B$ do not interact, and we want to see species $A$ and $B$ both
co-occurring and not interacting 10 times to be confident this is a true
negative, then we need an expected 1000 observations of all species if the
relative abundances of $A$ and $B$ are both $0.1$. 

Because the true FNR is latent, we can never actually be sure what the actual
number of false negatives in our data---however, we can use simulation to
estimate it for datasets of a given size using neutral models of observation. If
some of the “worst-case” FNRs presented in @fig:geometric(A) seem unrealistically
high, consider that species are observed in proportion to their relative
abundance. In the next section we demonstrate that the distribution of abundance
in ecosystems can lead to very high realized values of FNR ($p_{fn}$) simply as an
artifact of sampling effort.

![(A) The probability that an observed interaction is a true negative (y-axis)
given how many times it has been sampled as a non-interaction (x-axis). Each
color reflects a different value of $p_{fn}$, the false-negative rate
(FNR)---this is effectively the cdf of the geometric distribution. (B): The
expected needed observations of all individuals of all species (y-axis) required
to obtain a goal number of observations (colors) of a particular species, and a
function of the relative abundance of that focal species (x-axis). (C): False
negative rate (y-axis) as a function of total sampling effort (x-axis) and
network size, computed using the method described above. For 500 independent
draws from the niche model (@Williams2000SimRul) at varying levels of species
richness (colors) with connectance drawn according to the flexible-links model
(@MacDonald2020RevLin) as described in the main text. For each draw from the
niche model, 200 sets of 1500 observations are simulated, for which each the
mean false negative rate at each observation-step is
computed. Means denoted with points, with 1 in the first shade and 2 in the
second. (D): Same as (C), except using empirical food webs from Mangal database,
where richness. The outlier on (D) is a 714 species
food-web.](./figures/fig1.png){#fig:geometric}


# False-negatives as a product of relative abundance

We now show that the realized FNR changes drastically with sampling effort due
to the intrinsic variation of the abundance of individuals of each species
within a community. We do this by simulating the process of observation of
species interactions, applied both to 243 empirical food webs from the Mangal
database [@Banville2021ManJl] and random food-webs generated using the
niche model [@Williams2000SimRul], a simple generative model of food-web
structure. Our neutral model of observation assumes each observed species is
drawn from the distribution in proportion to each species' abundance at that
place and time.  The abundance distribution of a community can be
reasonably-well described by a log-normal distribution [@Volkov2003NeuThe]. In
addition to the log-normal distribution, we also tested the case where the
abundance distribution is derived from power-law scaling $Z^{(log(T_i)-1)}$
where $T_i$ is the trophic level of species $i$ and $Z$ is a scaling coefficient
[@Savage2004EffBod], which yields the same qualitative behavior. The practical
consequence of abundance distributions spanning many orders of magnitude of
abundance is that observing  two "rare" species interacting requires two low
probability events: observing two rare species at the same time.

To simulate the process of observation, for an ecological network $M$ with $S$
species, we sample abundances for each species from a standard-log-normal
distribution. For each true interaction in the adjacency matrix $M$ (i.e.
$M_{ij}=1$) we estimate the probability of observing both species $i$ and $j$ at
a given place and time by simulating $n$ observations of all individuals of any
species, where the species of the individual observed at the
$\{1,2,\dots,n\}$-th observation is drawn from the generated log-normal
distribution of abundances. For each pair of species $(i,j)$, if both $i$ and
$j$ are observed within the n-observations, the interaction is tallied as a true
positive if $M_{ij}=1$. If only one of $i$ or $j$ are observed---but not
both---in these $n$ observations, but $M_{ij}=1$, this is counted as a
false-negative, and a true-negative otherwise.


In @fig:geometric(C) we see this model of observation applied to niche model
networks across varying levels of species richness, and in @fig:geometric(D) the
observation model applied to Mangal food webs. For all niche model simulations
in this manuscript, for a given number of species $S$ the number of interactions
is drawn from the flexible-links model fit to Mangal data
[@MacDonald2020RevLin], effectively drawing the number of interactions $L$ for a
random niche model food-web as

$$L \sim  \text{BetaBinomial}(S^2-S+1,\mu\phi, 1-\mu\phi)$$

where the MAP estimate of $(\mu, \phi)$ applied to Mangal data from
[@MacDonald2020RevLin] is $(\mu=0.086, \phi=24.3)$. All simulations were done
with 500 independent replicates per unique number of observations $n$. All
analyses presented here are done in Julia v1.8 [@Bezanson2015JulFre] using both
EcologicalNetworks.jl v0.5 and Mangal.jl v0.4 [@Banville2021ManJl] and are
hosted at (GITHUB_LINK_TODO). Note that the empirical data, for the reasons
described above, very likely already contains many false negatives, we'll
revisit this issue in the final section.

From @fig:geometric(C) it is evident that the number of species considered in a
study is inseparable from the false-negative rate in that study, and this effect
should be taken into account when designing samples of ecological networks in
the future. We see a similar qualitative pattern in  @fig:geometric(D) where the
FNR drops off quickly as a function of observation effort, mediated by total
richness. The practical consequence of the bottom row of Figure 1 is whether the
total number of observations of all species (the x-axis) for the range of
possible FNR we deem acceptable (the y-axis) is feasible. This raises two
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
across the species pool? This argument is well-considered when sampling
individual species [@Willott2001SpeAcc], but has not yet been adopted for
designing samples of communities.

# Positive associations increase the false-negative rate

This model above doesn't consider the possibility that there are positive or
negative associations which shift the probability of observing two species
together due to their interaction [@Cazelles2016TheSpe]. However, here we
demonstrate that the probability of observing a false negative can be higher if
there is some positive association in the occurrence of species $A$ and $B$. If
we denote the probability that we observe the co-occurrence of two species $A$
and $B$ that we know interact as $P(AB)$ and if there is no association between
the marginal probabilities of observing $A$ and observing $B$, denoted $P(A)$
and $P(B)$ respectively, then the probability of observing their co-occurrence
$P(AB) = P(A)P(B)$. In the other case where there is some positive strength of
association between observing both $A$ and $B$ because this interaction is
"important" for each species, then the probability of observation both $A$ and
$B$, $P(AB)$, is greater than $P(A)P(B)$ as $P(A)$ and $P(B)$ are not
independent and instead are positively correlated, i.e. $P(AB)> P(A)P(B)$. In
this case, the probability of observing a false negative in our naive model from
@fig:geometric(A) is $p_{fn}= 1-P(AB)$, which due to the above inequality
implies $p_{fn} >1-P(A)P(B)$. This indicates an increasingly greater probability
of a false negative as $P(AB) \to P(AB) \gg P(A)P(B)$. However, this still does
not consider variation in species abundance in space and time
[@Poisot2015SpeWhy]. If positive or negative associations between species
structure variation in the distribution of $P(AB)$ across space/time, then the
spatial/temporal biases induced by data collection would further impact the
realized false negative rate, as the probability of false negative would not be
constant for each pair of species across sites. 

To test for this association in data we scoured Mangal for datasets with many
spatial or temporal replicates of the same system. For each dataset, we compute
the marginal probability $P(A)$ of occurrence of each species $A$ across all
networks in the dataset. For each pair of interacting species $A$ and $B$, we then
compute and compare the probability of co-occurrence if each species occurs
independently, $P(A)P(B)$, to the empirical joint probability of co-occurrence,
$P(AB)$. Following our analysis above, if $P(AB)$ is greater than $P(A)P(B)$, then we
expect our neutral estimates of the FNR above to underestimate the realized FNR.
In @fig:mangal, we see the difference between $P(AB)$ and $P(A)P(B)$ for the seven
suitable datasets with enough spatio-temporal replicates and a shared taxonomic
backbone (meaning all individual networks use common species identifiers) found
on Mangal to perform this analysis. Further details about each dataset are
reported in @tbl:id.  

![The difference between joint-probability of co-occurrence ($P(AB)$) and expected probability of co-occurrence under independence ($P(A)P(B)$) for interacting species for each dataset. The red-dashed line indicates 0. Each histogram represents a density, meaning the area of the entire curve sums to 1. The continuous density estimate (computed using local smoothing) is shown in grey. The p-value on each plot is the result of a one-sided t-test comparing the mean of each distribution to 0.](./figures/fig2.png){#fig:mangal}

|           Network          |      Type     |  $N$ |  $S$  |   $C$   |   $\bar{S}$   |   $\bar{C}$   |   $\bar{\beta}_{OS}$  |   $\bar{\beta}_{WN}$  |
|:--------------------------:|:-------------:|:---:|:---:|:-----:|:-----:|:-----:|:-----:|:-----:|
| @Kopelke2017FooStr   | Food Web      | 100 |  98 | 0.037 |  7.87 | 0.142 | 1.383 | 1.972 |
| @Thompson2000ResSol | Food Web      |  18 | 566 | 0.014 | 80.67 | 0.049 | 1.617 | 1.594 |
| @Havens1992ScaStr            | Food Web      |  50 | 188 | 0.065 | 33.58 | 0.099 | 1.468 | 1.881 |
| @Ponisio2017OppAtt      | Pollinator    | 100 | 226 | 0.079 |  23.0 | 0.056 | 1.436 | 1.870 |
| @Hadfield2014TalTwo     | Host-Parasite |  51 | 327 | 0.085 | 32.71 | 0.337 | 1.477 | 1.952 |
| @Closs1994SpaTem        | Food Web      |  12 |  61 |  0.14 | 29.09 | 0.080 | 1.736 | 1.864 |
| @CaraDonna2017IntRew    | Pollinator    |  86 | 122 |  0.18 | 21.42 | 0.312 | 1.527 | 1.907 |
Table: This table describes the datasets used in the above analysis (Fig 2). The table reports the type of each dataset, the total number of networks in each dataset $(N)$, the total species richness in each dataset $(S)$, the connectance of each metaweb (all interactions across the entire spatial-temporal extent) $(C)$, the mean species richness across each local network $S$, the mean connectance of each local network $C$, the mean $\beta$-diversity among overlapping species across all pairs of network species ($\bar{\beta}_{OS}$), and the mean $\beta$-diversity among all species in the metaweb ($\bar{\beta}_{WN}$). Both metrics are computed using KGL $\beta$-diversity [@Koleff2003MeaBet] {#tbl:id}


In each of these datasets, the joint probably of co-occurrence $P(AB)$ is
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

# The impact of false-negatives on network properties and prediction

Here, we assess the effect of false negatives on our ability to make predictions
about interactions, as well as their effect on network structure. The prevalence
of false-negatives in our data is the catalyst for interaction prediction in the
first place, and as a result methods have been proposed to counteract this bias 
[@Stock2017LinFil, @Poisot2022NetEmb]. However, it is feasible that the FNR
is so high that it could induce too much noise for an interaction prediction
model to detect the signal of possible interaction between species.

To test this we use the dataset from @Hadfield2014TalTwo that describes
host-parasite interaction networks sampled across 51 sites, and the same method
as @Strydom2021RoaPre to extract latent features for each species in this
dataset based on co-occurrence. We then predict a metaweb (equivalent to
predicting true or false for an interaction for each species pair, effectively a
binary classification problem) from these species-level features using four
candidate models for binary classification---three often used machine-learning
(ML) methods (Boosted Regression Tree (BRT), Random Forest (RF), Decision Tree
(DT)), and one naive model from classic statistics (Logistic Regression (LR)).
Each of the ML models are bootstrap aggregated (or bagged) with 100 replicates
each. We partition the data into 80-20 training-test split, and then seed the
training data with false negatives at varying rates, but crucially do nothing to
the test data. We fit all of these models using MLJ.jl, a high-level Julia
framework for a wide-variety of ML models [@Blaom2020MljJul]. We evaluate the
efficacy of these models using two common measures of binary classifier
performance: the area under the receiver-operator curve (ROC-AUC) and the area
under the precision-recall curve (PR-AUC), for more details see
[@Poisot2022GuiPre] Here, PR-AUC is slightly more relevant as it is a better
indicator of prediction of false-negatives. The results of these simulations are
shown in @fig:addedfnr(A) and (B).


![(A) The area-under the receiver-operator curve (ROC-AUC) and (B) The area-under the precision-recall curve (PR-AUC; right) for each different predictive model (colors/shapes) across a spectrum of the proportion of added false negatives (x-axis). (C) The mean trophic-level of all species in a network generated with the niche model across different species richnesses (colors). For each value of the FNR, the mean trophic level was computed across 50 replicates. The shaded region for each line is one standard-deviation across those replicates.](./figures/fig3.png){#fig:addedfnr}

One interesting result seen in @fig:addedfnr(A) and (B) is that the ROC-AUC
value does not approach random in the same way the PR-AUC curve does as we
increase the added FNR. The reason for this is that ROC-AUC is fundamentally not
as useful a metric in assessing predictive capacity as PR-AUC. As we keep adding
more false-negatives, the network eventually becomes a zeros matrix, and these
models can still learn to predict “no-interaction” for all possible species
pairs, which does far better than random guessing (ROC-AUC > 0.5) in terms of
the true and false positive rates (the components of ROC-AUC). This highlights a
more broad issue of label class imbalance, meaning there are far more
non-interactions than interactions in data. A full treatment of the importance
of class-balance is outside the scope of this paper, but is explored in-depth in
[@Poisot2022GuiPre].

Although these ML models are surprisingly performant at link prediction given
their simplicity, there have been several major developments in the field
applying deep-learning methods to many tasks in network inference and
prediction---namely graph-representation learning (GRL, @Khoshraftar2022SurGra)
and graph convolutional networks [@Zhang2019GraCon]. At this time, these
advances can not yet be applied to ecological networks because they require far
more data than we currently have. We already have lots of features that could be
used as inputs into these models (i.e. species level data about occurrence,
genomes, abundance, etc.), but our network datasets barely get into the hundreds
of local networks sampled across space and time @tbl:id. Once we start to get
into the thousands, these models will become more useful, but this can only be
done with systematic monitoring of interactions. This again highlights the need
to optimize our sampling effort to maximize the amount of information contained
in our data given the expensive nature of sampling interactions.  

We also consider how the FNR affects network properties. In @fig:addedfnr(C) we see
the mean trophic level across networks simulated using the niche model (as
above), across a spectrum of FNR values. In addition to the clear dependence on
richness, we see that mean trophic level, despite varying widely between niche
model simulations, tends to be relatively robust to false negatives and does not
deviate widely from the true value until very large FNRs, i.e. $p_{fn} > 0.7$. 

# Discussion 

Species interactions enable the persistence and functioning of ecosystems, but
our understanding of interactions is limited due to the intrinsic difficulty of
sampling. Here we have provided a null model for the expected number of
false-negatives in an interaction dataset. We demonstrated that we expect many
false-negatives in species interaction datasets purely due to the intrinsic
variation of abundances within a community. We also, for the first time to our
knowledge, measured the strength of association between co-occurrence and
interactions [@Cazelles2016TheSpe] across many empirical systems for the first
time, and found that these positive associations are both very common, and
showed algebraically that they increase the realized FNR.  We have also shown
that false-negatives could further impact our ability to both predict
interactions and infer properties of the networks, which highlights the need for
further research into methods for correcting this bias in existing data.

A better understanding of how false-negatives impact our inference of network
structure and dynamics is prediction of ecological networks by both using
information about individual species interactions and the structure of metawebs
on large scales is a practical necessity. False-negatives could pose a problem
for many forms of inference in network ecology. For example, inferring the
structural or dynamic stability of a network could be prone to error if the
observed network is not sampled "enough". What exactly "enough" means is then
specific to the application, and should be assessed via methods like those here
when designing samples. Further, predictions about network rewiring [@Thompson2017DisGov] due to range shifts in response to climate change could be
error-prone without accounting for interactions that have not been observed but
that still may become climatically infeasible. As is evident from Figure 1(A),
we can never guarantee there are no false-negatives in data. In recent years,
there has been interest toward explicitly accounting for false-negatives in
models [@Stock2017LinFil; @Young2021RecPla], and a predictive approach to
networks---rather than expecting our samples to fully capture all interactions
[@Strydom2021RoaPre]. As a result, better models for predicting interactions
are needed for interaction networks. This includes explicitly accounting for
observation error [@Johnson2021BayEst]---certain classes of models have
been used to reflect hidden states which account for detection error in
occupancy modeling [@Joseph2020NeuHie], and could be integrated in the predictive
models of interactions in the future.

A brief caveat here is that we do not consider the rate of false-positives---in
large part false-positives can be explained by misidentification of species,
although this could be a relevant consideration in some cases. The same logic
that we apply to false-negatives could easily be applied to false-positives,
e.g. that we can be much more confident that an interaction is a true positive
if we have observed it 50 times rather than only once, and we could similarly
model this using the geometric distribution as in @fig:geometric(A). However, because
ecological networks are so sparse, there are far more negatives than positives
in the dataset, and therefore likely to be far more false-negatives than
false-positives in absolute terms.  

This work has several practical consequences for the design of interaction
samples. Simulating the process of observation could be a powerful tool for
estimating the sampling effort required by a study that takes relative abundance
into account, and provides a null baseline for expected FNR. It is necessary to
take the size of the species pool into account when deciding how many total
samples is sufficient for an “acceptable” FNR (@fig:geometric(C & D)). Further the spatial
and temporal turnover of interactions means any approach to sampling
prioritization must be spatiotemporal. We demonstrated earlier that observed
negatives outside of the range of both species aren’t informative, and therefore
using species distribution models could aid in this spatial prioritization of
sampling sites.

Our work highlights the need for a quantitatively robust approach to sampling
design, both for interactions [@Jordano2016SamNet] and all other aspects of
biodiversity [@Carlson2020WhaWou]. As anthropogenic forces create rapid shifts
in our planet's climate and biosphere, this is an imperative to maximize the
amount of ecological information we get in our finite samples, and make our
inferences and decisions based on this data as robust as possible. Where we
choose to sample, and how often we choose to sample there, has strong impacts on
the inferences we make from data. Incorporating a better understanding of
sampling effort and bias to the design of biodiversity monitoring systems, and
the inference and predictive models we apply to this data, is imperative in
understanding how biodiversity is changing, and making actionable forecasts
about the future of ecological interactions on our planet.


# References
