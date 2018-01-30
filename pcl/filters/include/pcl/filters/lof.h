/*
* Software License Agreement (BSD License)
*
*  Point Cloud Library (PCL) - www.pointclouds.org
*  Copyright (c) 2018-2020, Francesco Pirotti, CIRGEO University of Padova, <francesco.pirotti@unipd.it>.
*
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of Willow Garage, Inc. nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*/


#ifndef PCL_FILTERS_LOF_H_
#define PCL_FILTERS_LOF_H_

#include <pcl/filters/filter_indices.h>
#include <pcl/search/pcl_search.h> 

namespace pcl
{
	/** \brief @b LOFFilter uses point neighborhood statistics to filter outlier data.
	* \details The algorithm iterates through the entire input twice:
	* During the next iteration the points will be classified as inlier or outlier if their average neighbor distance is below or above this threshold respectively.
	* <br>
	* The neighbors found for each query point will be found amongst ALL points of setInputCloud(), not just those indexed by setIndices().
	* The setIndices() method only indexes the points that will be iterated through as search query points.
	* <br><br>
	* For more information:  Breunig MM, Kriegel H-P, Ng RT, Sander J. 
	* LOF: Identifying Density-Based Local Outliers. Proc. 2000 Acm Sigmod Int. Conf. Manag. Data [Internet]. 2000;1–12. 
	* Available from: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.35.8948
	* <br><br>
	* Usage example:
	* \code
	* pcl::LOFFilter<PointType> sorfilter (true); // Initializing with true will allow us to extract the removed indices
	* sorfilter.setInputCloud (cloud_in);
	* sorfilter.setK (20);
	* sorfilter.setThresh (1.0);
	* sorfilter.filter (*cloud_out);
	* // The resulting cloud_out contains all points of cloud_in that have an average distance to their 8 nearest neighbors that is below the computed threshold
	* // Using a standard deviation multiplier of 1.0 and assuming the average distances are normally distributed there is a 84.1% chance that a point will be an inlier
	* indices_rem = sorfilter.getRemovedIndices ();
	* // The indices_rem array indexes all points of cloud_in that are outliers
	* \endcode
	* \author Francesco Pirotti
	* \ingroup filters
	*/
	struct pt
	{
		//float kDistance;
		float lrd;
	//	float lof;
		std::vector<int> knn;
		std::vector<float> kdists;
	};
	template<typename PointT>
	class LOFFilter : public FilterIndices<PointT>
	{
	protected:
		typedef typename FilterIndices<PointT>::PointCloud PointCloud;
		typedef typename PointCloud::Ptr PointCloudPtr;
		typedef typename PointCloud::ConstPtr PointCloudConstPtr;
		typedef typename pcl::search::Search<PointT>::Ptr SearcherPtr;

	public:

		typedef boost::shared_ptr< LOFFilter<PointT> > Ptr;
		typedef boost::shared_ptr< const LOFFilter<PointT> > ConstPtr;


		/** \brief Constructor.
		* \param[in] extract_removed_indices Set to true if you want to be able to extract the indices of points being removed (default = false).
		*/
		LOFFilter(bool extract_removed_indices = false) :
			FilterIndices<PointT>::FilterIndices(extract_removed_indices),
			searcher_(),
			k_(20), 
			threshold_(2.0),
			lof_(new std::vector<float>)
		{
			filter_name_ = "LOFFilter";
		}

		/** \brief Set the number of nearest neighbors to use for local reachability distance estimation.
		* \param[in] nr_k The number of points to use for mean distance estimation.
		*/
		inline void
			setK(int nr_k)
		{
			k_ = nr_k;
		}

		/** \brief Get the number of nearest neighbors to use for local reachability distance estimation.
		* \return The number of points to use for mean distance estimation.
		*/
		inline int
			getK()
		{
			return (k_);
		}

		/** \brief Set threshold of LOF.
		* \details .
		* Points will be classified as inlier or outlier if their average neighbor distance is below or above this threshold respectively.
		* \param[in] threshold_ The threshold.
		*/
		inline void
			setThresh(double threshold)
		{
			threshold_ = threshold;
		}

		/** \brief Get the threshold of LOF.
		* \details The LOF value threshold will be equal to: threshold .
		* Points will be classified as inlier or outlier if their LOF Score is below or above this threshold respectively.
		*/
		inline double
			getThresh()
		{
			return (threshold_);
		}


		/** \brief Get the vector of LOF scores.
		* \details The vector<float> of LOF values for each  point.
		*/
		inline boost::shared_ptr <std::vector<float> >
			getLOFs()
		{
			return (lof_);
		}
	protected:
		using PCLBase<PointT>::input_; 
		using PCLBase<PointT>::indices_;
		using Filter<PointT>::filter_name_;
		using Filter<PointT>::getClassName;
		using FilterIndices<PointT>::negative_;
		using FilterIndices<PointT>::keep_organized_;
		using FilterIndices<PointT>::user_filter_value_;
		using FilterIndices<PointT>::extract_removed_indices_;
		using FilterIndices<PointT>::removed_indices_;

		/** \brief Filtered results are stored in a separate point cloud.
		* \param[out] output The resultant point cloud.
		*/
		void
			applyFilter(PointCloud &output);

		/** \brief Filtered results are indexed by an indices array.
		* \param[out] indices The resultant indices.
		*/
		void
			applyFilter(std::vector<int> &indices)
		{
			applyFilterIndices(indices);
		}

		/** \brief Filtered results are indexed by an indices array.
		* \param[out] indices The resultant indices.
		*/
		void
			applyFilterIndices(std::vector<int> &indices);

	private:
		/** \brief A pointer to the spatial search object. */
		SearcherPtr searcher_;
		 
		/** \brief The number of points to use for k-distance calculation. */
		int k_;

		/** \brief Pointer to vector of lof scores. */
		boost::shared_ptr <std::vector<float> > lof_;
		/** \brief Array of k-distances. */
		std::vector<float> kdistance_;
		/** \brief LOF Score threshold (i.e., points outside 
		* will be marked as outliers). */
		double threshold_;
	};

	/** \brief @b LOFFilter uses point neighborhood density statistics to filter outlier data. For more
	* information check:
	* Breunig, M. M., Kriegel, H.-P., Ng, R. T., & Sander, J. (2000).
	* LOF: Identifying Density-Based Local Outliers. Proceedings of the 2000 Acm Sigmod International Conference on Management of Data, 1–12.
	* http://doi.org/10.1145/335191.335388
	* 
	* \author Francesco Pirotti
	* \ingroup filters
	*/
	template<>
	class PCL_EXPORTS LOFFilter<pcl::PCLPointCloud2> : public FilterIndices<pcl::PCLPointCloud2>
	{
		using FilterIndices<pcl::PCLPointCloud2>::filter_name_;
		using FilterIndices<pcl::PCLPointCloud2>::getClassName;

		using FilterIndices<pcl::PCLPointCloud2>::removed_indices_;
		using FilterIndices<pcl::PCLPointCloud2>::extract_removed_indices_;

		typedef pcl::search::Search<pcl::PointXYZ> KdTree;
		typedef pcl::search::Search<pcl::PointXYZ>::Ptr KdTreePtr;

		typedef pcl::PCLPointCloud2 PCLPointCloud2;
		typedef PCLPointCloud2::Ptr PCLPointCloud2Ptr;
		typedef PCLPointCloud2::ConstPtr PCLPointCloud2ConstPtr;

	public:
		/** \brief Empty constructor. */
		LOFFilter(bool extract_removed_indices = false) :
			FilterIndices<pcl::PCLPointCloud2>::FilterIndices(extract_removed_indices), 
			k_(20),
			threshold_(2.0), 
			tree_(),
			lof_(new std::vector<float>)
		{
			filter_name_ = "LOFFilter";
		}

		/** \brief Set the number of points (k) to use for mean distance estimation
		* \param nr_k the number of points to use for mean distance estimation
		*/
		inline void
			setK(int nr_k)
		{
			k_ = nr_k;
		}

		/** \brief Get the number of points to use for mean distance estimation. */
		inline int
			getK()
		{
			return (k_);
		}


		/** \brief Set the standard deviation multiplier threshold. All points outside the
		* \f[ \mu \pm \sigma \cdot std\_mul \f]
		* will be considered outliers, where \f$ \mu \f$ is the estimated mean,
		* and \f$ \sigma \f$ is the standard deviation.
		* \param threshold the standard deviation multiplier threshold
		*/
		inline void
			setThresh(double threshold)
		{
			threshold_ = threshold;
		}

		/** \brief Get the standard deviation multiplier threshold as set by the user. */
		inline double
			getThresh()
		{
			return (threshold_);
		}

	protected:
		/** \brief The number of points to use for LOF calculation. */
		int k_;


		/** \brief Shared pointer to vector of lof scores. */ 
		boost::shared_ptr <std::vector<float> > lof_;
		/** \brief Array of k-distances. */
		std::vector<float> kdistance_;

		/** \brief Standard deviations threshold (i.e., points outside of
		* \f$ \mu \pm \sigma \cdot std\_mul \f$ will be marked as outliers).
		*/
		double threshold_;

		/** \brief A pointer to the spatial search object. */
		KdTreePtr tree_;

		virtual void
			applyFilter(std::vector<int> &indices);

		virtual void
			applyFilter(PCLPointCloud2 &output);

		/**
		* \brief Compute the statistical values used in both applyFilter methods.
		*
		* This method tries to avoid duplicate code.
		*/
		virtual void
			generateStatistics();
	};
}

#ifdef PCL_NO_PRECOMPILE
#include <pcl/filters/impl/lof.hpp>
#endif

#endif // PCL_FILTERS_LOF_H_