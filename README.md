# Analysis of Drugnet to reduce drug spread and understand local structures in the network

In this project “Drugnet” network has been analyzed to find important nodes as well as underlaying structures in the network. I found 
that the network is male dominated and high percentage of participants belongs from Latino ethnicity. To find popularity of participants 
in the network in-degree centrality is used. From in-degree centrality I observed that there are many nodes which are isolates and
then decided to focus on nodes which have in-degree centrality between 10 to 5. So, there are total 12 nodes which I selected and if we 
educate them about drug usage then there is large probability that more than 50% of the network can get influenced. Further, I also used
ERGMs for understanding the underlaying structures in the network and compared simulated models with observed network. From this analysis
I found that edges, mutuality, asymmetric dyads, isolates, homophily and simmelian triads are the important local structures which 
represents some part of the global pattern in the network. 

# References
1)	WEEKS, M. R., CLAIR, S., BORGATTI, S. P., RADDA, K. & SCHENSUL, J. J. 2002. Social networks of drug users in high-risk sites: Finding the connections. AIDS and Behaviour, 6, 193-206.
2)	Introduction to Exponential-family Random Graph (ERG or p*) modeling with ergm, The Statnet Devlopment Team, August 2017
3)	ergm: A Package to Fit, Simulate and Diagnose Exponential-Family Models for Networks, Journal of Statistical Software, David R. Hunter et. al. (2008)
4)	Speciation of Exponential-Family Random Graph Models: Terms and Computational Aspects, Journal of Statistical Software, Martina Morris et. al. (2008)
5)	Social Network Analysis with sna, Journal of Statistical Software, Carter T. Butts (2008)
6)	Exponential Random Graph (ERG or p*) Models, COMM 645: Communication Networks Annenberg School of Communication University of Southern California
7)	A statnet Tutorial, Journal of Statistical Software, Steven M. Goodreau et. al. (2008)
8)	https://rstudio-pubs-static.s3.amazonaws.com/157501_93a72a58ec614946901e10edf78c1384.html
9)	https://www.sci.unich.it/~francesc/teaching/network/
10)	http://badhessian.org/2012/09/lessons-on-exponential-random-graph-modeling-from-greys-anatomy-hook-ups/
11)	http://www.mjdenny.com/Preparing_Network_Data_In_R.html
