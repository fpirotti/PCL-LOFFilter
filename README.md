# PCL-LOFFilter
For more information:  Breunig MM, Kriegel H-P, Ng RT, Sander J. 
	 LOF: Identifying Density-Based Local Outliers. Proc. 2000 Acm Sigmod Int. Conf. Manag. Data [Internet]. 2000;1–12. 
	 Available from: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.35.8948
	* Usage example:
	 pcl::LOFFilter<PointType> sorfilter (true); // Initializing with true will allow us to extract the removed indices
	 sorfilter.setInputCloud (cloud_in);
	 sorfilter.setK (20);
	 sorfilter.setThresh (1.0);
	 sorfilter.filter (*cloud_out);
