/*
* Software License Agreement (BSD License)
*
*  Point Cloud Library (PCL) - www.pointclouds.org
*  Copyright (c) 2018-, Francesco Pirotti, CIRGEO University of Padova, <francesco.pirotti@unipd.it>.
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
*   * Neither the name of the copyright holder(s) nor the names of its
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
* $Id$
*
*/

#include <pcl/filters/impl/lof.hpp>
#include <pcl/conversions.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////
void
pcl::LOFFilter<pcl::PCLPointCloud2>::applyFilter(PCLPointCloud2 &output)
{

	PCL_INFO("[pcl::%s::applyFilter] 22 Input/Output dataset do not have scalar field \"strength\" from pcl::InterestPoint - please create an input point cloud of type pcl::InterestPoint so that LOF value can be saved in the scalar field !\n", getClassName().c_str());
	PCL_WARN("[pcl::%s::applyFilter] 22 Input/Output dataset do not have scalar field \"strength\" from pcl::InterestPoint - please create an input point cloud of type pcl::InterestPoint so that LOF value can be saved in the scalar field !\n", getClassName().c_str());

	bool scalar_not_available = getFieldIndex(output, "strength") < 0;
	// If field  strength is not present
	if ( scalar_not_available)
	{
		PCL_ERROR("[pcl::%s::applyFilter] Input/Output dataset do not have scalar field \"strength\" from pcl::InterestPoint - please create an input point cloud of type pcl::InterestPoint so that LOF value can be saved in the scalar field !\n", getClassName().c_str());
		output.width = output.height = 0;
		output.data.clear();
		return;
	}

	// If fields x/y/z are not present, we cannot filter
	if (x_idx_ == -1 || y_idx_ == -1 || z_idx_ == -1)
	{
		PCL_ERROR("[pcl::%s::applyFilter] Input dataset doesn't have x-y-z coordinates!\n", getClassName().c_str());
		output.width = output.height = 0; 
		output.data.clear();
		return;
	}

	if (threshold_ == 0.0)
	{
		PCL_ERROR("[pcl::%s::applyFilter] LOF Score threshold %.2f not correct!\n", getClassName().c_str(), threshold_);
		output.width = output.height = 0;
		output.data.clear();
		return;
	}
	 
	generateStatistics();
//	double const threshold_ = mean + threshold_ * stddev; // a distance that is bigger than this signals an outlier

																// Copy the common fields
	output.is_dense = input_->is_dense;
	output.is_bigendian = input_->is_bigendian;
	output.point_step = input_->point_step;
	if (keep_organized_)
	{
		output.width = input_->width;
		output.height = input_->height;
		output.data.resize(input_->data.size());
	}
	else
	{
		output.height = 1;
		output.data.resize(indices_->size() * input_->point_step); // reserve enough space
	}

	removed_indices_->resize(input_->data.size());


	// Build a new cloud by neglecting outliers
	int nr_p = 0;
	int nr_removed_p = 0;
	bool remove_point = false;
	for (int cp = 0; cp < static_cast<int> (indices_->size()); ++cp)
	{


		if (negative_)
			remove_point = ((*lof_)[cp] <= threshold_);
		else
			remove_point = ((*lof_)[cp] > threshold_);
		
		if (remove_point)
		{
			if (extract_removed_indices_)
				(*removed_indices_)[nr_removed_p++] = cp;

			if (keep_organized_)
			{
				/* Set the current point to NaN. */
				*(reinterpret_cast<float*>(&output.data[nr_p * output.point_step]) + 0) = std::numeric_limits<float>::quiet_NaN();
				*(reinterpret_cast<float*>(&output.data[nr_p * output.point_step]) + 1) = std::numeric_limits<float>::quiet_NaN();
				*(reinterpret_cast<float*>(&output.data[nr_p * output.point_step]) + 2) = std::numeric_limits<float>::quiet_NaN();
				nr_p++;
				output.is_dense = false;
			}
			else
				continue;
		}
		else
		{
			memcpy(&output.data[nr_p * output.point_step], &input_->data[(*indices_)[cp] * output.point_step],
				output.point_step);
			nr_p++;
		}
	}

	if (!keep_organized_)
	{
		output.width = nr_p;
		output.data.resize(output.width * output.point_step);
	}
	output.row_step = output.point_step * output.width;

	removed_indices_->resize(nr_removed_p);
}

///////////////////////////////////////////////////////////////////////////////////////////
void
pcl::LOFFilter<pcl::PCLPointCloud2>::applyFilter(vector<int>& indices)
{
	PCL_INFO("[pcl::%s::applyFilter] 33 Input/Output dataset do not have scalar field \"strength\" from pcl::InterestPoint - please create an input point cloud of type pcl::InterestPoint so that LOF value can be saved in the scalar field !\n", getClassName().c_str());
	PCL_WARN("[pcl::%s::applyFilter] 223 Input/Output dataset do not have scalar field \"strength\" from pcl::InterestPoint - please create an input point cloud of type pcl::InterestPoint so that LOF value can be saved in the scalar field !\n", getClassName().c_str());

	// If fields x/y/z are not present, we cannot filter
	if (x_idx_ == -1 || y_idx_ == -1 || z_idx_ == -1)
	{
		PCL_ERROR("[pcl::%s::applyFilter] Input dataset doesn't have x-y-z coordinates!\n", getClassName().c_str());
		indices.clear();
		return;
	}

	if (threshold_ == 0.0)
	{
		PCL_ERROR("[pcl::%s::applyFilter] Standard deviation multipler not set!\n", getClassName().c_str());
		indices.clear();
		return;
	}

	generateStatistics();
	//double const threshold_ = mean + threshold_ * stddev; // a distance that is bigger than this signals an outlier

																// Second pass: Classify the points on the computed distance threshold
	size_t nr_p = 0, nr_removed_p = 0;
	for (size_t cp = 0; cp < indices_->size(); ++cp)
	{
		// Points having a too high average distance are outliers and are passed to removed indices
		// Unless negative was set, then it's the opposite condition
		if ((!negative_ && (*lof_)[cp] > threshold_) || (negative_ && (*lof_)[cp] <= threshold_))
		{
			if (extract_removed_indices_)
				(*removed_indices_)[nr_removed_p++] = (*indices_)[cp];
			continue;
		}

		// Otherwise it was a normal point for output (inlier)
		indices[nr_p++] = (*indices_)[cp];
	}

	// Resize the output arrays
	indices.resize(nr_p);
	removed_indices_->resize(nr_p);
}

///////////////////////////////////////////////////////////////////////////////////////////
void
pcl::LOFFilter<pcl::PCLPointCloud2>::generateStatistics()
{
	// Send the input dataset to the spatial locator
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::fromPCLPointCloud2(*input_, *cloud);

	pt *ptn = new pt[indices_->size()];

	// Initialize the spatial locator
	if (!tree_)
	{
		if (cloud->isOrganized())
			tree_.reset(new pcl::search::OrganizedNeighbor<pcl::PointXYZ>());
		else
			tree_.reset(new pcl::search::KdTree<pcl::PointXYZ>(true));
	}

	tree_->setInputCloud(cloud);
	// Allocate enough space to hold the results  
 
	(*lof_).resize(indices_->size());
	 
	int valid_distances = 0;
	// TWO loops
	// Go over all the points and calculate the k-distance
	for (size_t cp = 0; cp < indices_->size(); ++cp)
	{

		if (!pcl_isfinite(cloud->points[(*indices_)[cp]].x) ||
			!pcl_isfinite(cloud->points[(*indices_)[cp]].y) ||
			!pcl_isfinite(cloud->points[(*indices_)[cp]].z))
		{
			(*lof_)[cp] = 0;
			continue;
		}
		 
		// Minimum distance (if k_ == 2) or mean distance

		ptn[cp].knn.resize(k_);
		ptn[cp].kdists.resize(k_);
		if (tree_->nearestKSearch((*indices_)[cp], k_, ptn[cp].knn, ptn[cp].kdists) == 0)
		{
			(*lof_)[cp] = 0;
			PCL_WARN("[pcl::%s::applyFilter] Searching for the closest %d neighbors failed.\n", getClassName().c_str(), k_);
			continue;
		}

	}


	for (size_t cp = 0; cp < indices_->size(); ++cp)
	{
		double summedReachDist = 0;

		//Iterate over the K neighbors to ith point
		for (int j = 0; j < k_; j++)
		{
			int k = ptn[cp].knn[j];
			pt p2 = ptn[ k ];
			// dists are sorted so last is the farther one
			//float kDist = p2.kdists[k_-1];
			float reachDist = max(p2.kdists[(k_ - 1)], ptn[cp].kdists[j] );

			summedReachDist += reachDist;
		}

		ptn[cp].lrd = (float)k_ / summedReachDist;
	}

	for (size_t cp = 0; cp < indices_->size(); ++cp)
	{
		float summedLRDRatio = 0;

		//Iterate over the MINPTS neighbors to ith point
		for (int j = 0; j < k_; j++)
		{
			int k = ptn[cp].knn[j];
			pt p2 = ptn[k];
			// dists are sorted so last is the farther one
			//float kDist = p2.kdists[k_-1];
			float reachDist = max(p2.kdists[(k_ - 1)], ptn[cp].kdists[j]);

			summedLRDRatio += (p2.lrd / ptn[cp].lrd);
			 
		}

		(*lof_)[cp] = (float)summedLRDRatio / (float)k_;
	}

}


#ifndef PCL_NO_PRECOMPILE
#include <pcl/impl/instantiate.hpp>
#include <pcl/point_types.h>

// Instantiations of specific point types
PCL_INSTANTIATE(LOFFilter, PCL_XYZ_POINT_TYPES)

#endif    // PCL_NO_PRECOMPILE

