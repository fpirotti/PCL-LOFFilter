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



#ifndef PCL_FILTERS_LOF_IMPL_H_
#define PCL_FILTERS_LOF_IMPL_H_

#include <pcl/filters/lof.h>

#include <pcl/common/io.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT> void
pcl::LOFFilter<PointT>::applyFilter(PointCloud &output)
{

	std::vector<int> indices;
	if (keep_organized_)
	{
		bool temp = extract_removed_indices_;
		extract_removed_indices_ = true;
		applyFilterIndices(indices);
		extract_removed_indices_ = temp;

		output = *input_;
		for (int rii = 0; rii < static_cast<int> (removed_indices_->size()); ++rii)  // rii = removed indices iterator
			output.points[(*removed_indices_)[rii]].x = output.points[(*removed_indices_)[rii]].y = output.points[(*removed_indices_)[rii]].z = user_filter_value_;
		if (!pcl_isfinite(user_filter_value_))
			output.is_dense = false;
	}
	else
	{
		applyFilterIndices(indices);
		copyPointCloud(*input_, indices, output);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT> void
pcl::LOFFilter<PointT>::applyFilterIndices(std::vector<int> &indices)
{


	pt *ptn = new pt[indices_->size()];


	// Initialize the search class
	if (!searcher_)
	{
		if (input_->isOrganized())
			searcher_.reset(new pcl::search::OrganizedNeighbor<PointT>());
		else
			searcher_.reset(new pcl::search::KdTree<PointT>(true));
	}
	searcher_->setInputCloud(input_);

	// The arrays to be used
	std::vector<int> nn_indices(k_);
	std::vector<float> nn_dists(k_);
	//std::vector<float> distances(indices_->size());
	indices.resize(indices_->size());
	(*lof_).resize(indices_->size());
	removed_indices_->resize(indices_->size());
	int oii = 0, rii = 0;  // oii = output indices iterator, rii = removed indices iterator


	double summedReachDist = 0;
	PCL_INFO("[pcl::%s::applyFilter] First pass to get closest %d neighbors.\n", getClassName().c_str(), k_);
	for (int cp = 0; cp < static_cast<int> (indices_->size()); ++cp)  // cp = input indices iterator
	{
		if (!pcl_isfinite(input_->points[(*indices_)[cp]].x) ||
			!pcl_isfinite(input_->points[(*indices_)[cp]].y) ||
			!pcl_isfinite(input_->points[(*indices_)[cp]].z))
		{ 
			continue;
		}
		ptn[cp].knn.resize(k_);
		ptn[cp].kdists.resize(k_);
 
		// Perform the nearest k search
		if (searcher_->nearestKSearch((*indices_)[cp], k_ + 1, ptn[cp].knn, ptn[cp].kdists) == 0)
		{
			(*lof_)[cp] = 0;
			PCL_WARN("[pcl::%s::applyFilter] Searching for the closest %d neighbors failed.\n", getClassName().c_str(), k_);
			continue;
		}

	}

	// Second pass: Calculate local reachability distance

	PCL_INFO("[pcl::%s::applyFilter] Second pass to get Local Reachability Distance.\n", getClassName().c_str());
	for (int cp = 0; cp < static_cast<int> (indices_->size()); ++cp)  // cp = input indices iterator
	{


		double summedReachDist = 0;

		//Iterate over the K neighbors to ith point
		for (int j = 0; j < k_; j++)
		{
			int k = ptn[cp].knn[j];
			pt p2 = ptn[k];
			// dists are sorted so last is the farther one
			//float kDist = p2.kdists[k_-1];
			float reachDist = max(p2.kdists[(k_ - 1)], ptn[cp].kdists[j]);

			summedReachDist += reachDist;
		}

		ptn[cp].lrd = (float)k_ / summedReachDist;
	}

	PCL_INFO("[pcl::%s::applyFilter] Third pass to get Local Outlier Factor and apply threshold.\n", getClassName().c_str());
	// Third pass: Compute LOF and apply threshold
	for (int cp = 0; cp < static_cast<int> (indices_->size()); ++cp)  // cp = input indices iterator
	{

		float summedLRDRatio = 0;

		//Iterate over the MINPTS neighbors to ith point
		for (int j = 0; j < k_; j++)
		{
			int k = ptn[cp].knn[j];
			pt p2 = ptn[k];
			summedLRDRatio += (p2.lrd / ptn[cp].lrd);
		}

		(*lof_)[cp] = (float)summedLRDRatio / (float)k_;
		// Points having a too high average distance are outliers and are passed to removed indices
		// Unless negative was set, then it's the opposite condition
		if ((!negative_ && (*lof_)[cp]  > threshold_) || (negative_ && (*lof_)[cp]  <= threshold_))
		{
			if (extract_removed_indices_)
				(*removed_indices_)[rii++] = (*indices_)[cp];
			continue;
		}

		// Otherwise it was a normal point for output (inlier)
		indices[oii++] = (*indices_)[cp];
	}

	// Resize the output arrays
	indices.resize(oii);
	removed_indices_->resize(rii);
}


#define PCL_INSTANTIATE_LOFFilter(T) template class PCL_EXPORTS pcl::LOFFilter<T>;


#endif // PCL_FILTERS_LOF_H_