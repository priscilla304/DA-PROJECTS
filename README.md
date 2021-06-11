# Data Science Projects


**CORONA VIRUS CASES WORLD DATA - DATA VISUALIZATION:**

This project focuses on tracking the spread of coronavirus disease (COVID-19) and helps us visualize the areas affected by the virus around the globe. Tracking the COVID-19 spread enables us to understand the seriousness of the current situation. The global pandemic can be visualized using various types of plots/charts. The visualizations enable people to understand the data more quickly. We have used various techniques and color schemes to develop Coronavirus Dashboard to help visualize available data in an informative way.

DATA:

Data acquisition is the process of fetching the data from a source to our system. We used R for data acquisition from two different sources. The data for the Coronavirus Dashboard is taken from the below sources:
•	Covid19.analytics (R package)
https://cran.r-project.org/web/packages/covid19.analytics/vignettes/covid19.analytics.html
•	The New York Times (Github repository)
https://github.com/nytimes/covid-19-data/

**Twitter Sentiment Analysis - Data Mining2:**
Sentiment Analysis is one the interesting topic in Data Mining where we get the key words from the tweets and determine the public reaction. This helps not only helps us make decisions but also get a quick idea from a large set of data.
The steps followed in this projects are :
Extract Tweets using Keyword(Data Extraction Using Twitter Development Account)
Preprocessing the data using qdap package in R to remove url , hash-tag, emot-icon, email, digits, punctuations, stop words and special characters.
Creating bag of words and word-cloud.
Topic Modelling
Sentiment Analysis
Comparison cloud for the combined tweets of two set of tweets.
Text Clustering

**Particle Classification - Data Mining1:**
The goal in this task is to learn a classification rule that differentiates between two types of particles generated in high energy collider experiments. It is a binary classification problem with 78 attributes. The data set has 50,000 examples. Our task is to build a classification model that optimizes ROC area (maximize).

**CAPESTONE PROJECT:**
Image Segmentation using CT Scan data

The scope and the intent for the GlassJaw project for the Minimum Viable Product (MVP) is to use Machine Learning as a tool to locate and auto-segment medical scans to help assist with VSP’s and design implants for cranio-maxillofacial surgeries. Our team has implemented an unsupervised learning approach to create labels using Python by adopting methods such as k-means clustering. This method uses thresholding to filter out the mandible part of the DICOM (Digital Imaging and Communication in Medicine) to locate the specific area of interest such as a particular region of the bone through labeling each segmented part. The team had also gone further and developed on this by using automatic segmentation based supervised learning through convolutional neural networks (CNN).
