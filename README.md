# simple virtual screening
## A simple python script for ligand-based virtual screening

virtual_screening.py is a simple python script that I wrote to run ligand-based virtual screening calculations.

**The purpose of this script is limited to exercise python. For extensively tested tools, please, consult specialized scientific literature.**

Given a reference set and and a database of small molecule sencoded in fingerprints, the script computes a similarity (or distance) measure and returns the ranked database.

The script can compute rankings based on one of three different similarity measures: 

* Tanimoto index
* Dice index
* Cosine coefficient

or a distance measure

* Soergel distance

for a definition of these measures and a discussion of their utility, please, refer to reference 1.

#### References

1. Bajusz, D., Rácz, A. & Héberger, K. Why is Tanimoto index an appropriate choice for fingerprint-based similarity calculations? Journal of Cheminformatics 7, 20 (2015).
