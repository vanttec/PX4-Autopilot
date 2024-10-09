/****************************************************************************
 *
 *   Copyright (c) 2024 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

#pragma once

#include "mathlib/math/Limits.hpp"
#include <matrix/math.hpp>

class LatLonAlt
{
public:
	LatLonAlt() = default;
	LatLonAlt(const double latitude_deg, const double longitude_deg,
		  const float altitude_m) : _latitude(math::radians(latitude_deg)), _longitude(math::radians(longitude_deg)),
		_altitude(altitude_m) {};

	double latitude_deg() const { return math::degrees(_latitude); }
	double longitude_deg() const { return math::degrees(_longitude); }

	const double &latitude_rad() const { return _latitude; }
	const double &longitude_rad() const { return _longitude; }
	const float &altitude() const {return _altitude; }

	double &latitude_rad() { return _latitude; }
	double &longitude_rad() { return _longitude; }
	float &altitude() {return _altitude; }

	void operator+=(const matrix::Vector3f &delta_pos)
	{
		matrix::Vector2d d_lat_lon_to_d_xy = deltaLatLonToDeltaXY(_latitude, _altitude);
		_latitude = matrix::wrap_pi(_latitude + static_cast<double>(delta_pos(0)) / d_lat_lon_to_d_xy(0));
		_longitude = matrix::wrap_pi(_longitude + static_cast<double>(delta_pos(1)) / d_lat_lon_to_d_xy(1));
		_altitude -= delta_pos(2);
	}

	void operator+=(const matrix::Vector2f &delta_pos)
	{
		matrix::Vector2d d_lat_lon_to_d_xy = deltaLatLonToDeltaXY(_latitude, _altitude);
		_latitude = matrix::wrap_pi(_latitude + static_cast<double>(delta_pos(0)) / d_lat_lon_to_d_xy(0));
		_longitude = matrix::wrap_pi(_longitude + static_cast<double>(delta_pos(1)) / d_lat_lon_to_d_xy(1));
	}

	LatLonAlt operator+(const matrix::Vector3f &delta_pos) const
	{
		LatLonAlt lla_new{};
		matrix::Vector2d d_lat_lon_to_d_xy = deltaLatLonToDeltaXY(_latitude, _altitude);
		lla_new.latitude_rad() = matrix::wrap_pi(_latitude + static_cast<double>(delta_pos(0)) / d_lat_lon_to_d_xy(0));
		lla_new.longitude_rad() = matrix::wrap_pi(_longitude + static_cast<double>(delta_pos(1)) / d_lat_lon_to_d_xy(1));
		lla_new.altitude() = _altitude - delta_pos(2);
		return lla_new;
	}

	matrix::Vector3f operator-(const LatLonAlt &lla) const
	{
		const double delta_lat = matrix::wrap_pi(_latitude - lla.latitude_rad());
		const double delta_lon = matrix::wrap_pi(_longitude - lla.longitude_rad());
		const float delta_alt = _altitude - lla.altitude();

		const matrix::Vector2d d_lat_lon_to_d_xy = deltaLatLonToDeltaXY(_latitude, _altitude);
		return matrix::Vector3f(static_cast<float>(delta_lat * d_lat_lon_to_d_xy(0)),
					static_cast<float>(delta_lon * d_lat_lon_to_d_xy(1)),
					-delta_alt);
	}

private:
	// Convert between curvilinear and cartesian errors
	static matrix::Vector2d deltaLatLonToDeltaXY(const double latitude, const float altitude)
	{
		double r_n;
		double r_e;
		computeRadiiOfCurvature(latitude, r_n, r_e);
		const double dn_dlat = r_n + static_cast<double>(altitude);
		const double de_dlon = (r_e + static_cast<double>(altitude)) * cos(latitude);

		return matrix::Vector2d(dn_dlat, de_dlon);
	}

	static void computeRadiiOfCurvature(const double latitude, double &meridian_radius_of_curvature,
					    double &transverse_radius_of_curvature)
	{
		const double tmp = 1.0 - pow(Wgs84::eccentricity * sin(latitude), 2);
		const double sqrt_tmp = sqrt(tmp);
		meridian_radius_of_curvature = Wgs84::meridian_radius_of_curvature_numerator / (tmp * tmp * sqrt_tmp);
		transverse_radius_of_curvature = Wgs84::equatorial_radius / sqrt_tmp;
	}

	struct Wgs84 {
		static constexpr double equatorial_radius = 6378137;
		static constexpr double eccentricity = 0.0818191908425;
		static constexpr double meridian_radius_of_curvature_numerator = equatorial_radius * (1.0 - pow(eccentricity, 2));
	};

	double _latitude{NAN};  ///< (rad)
	double _longitude{NAN}; ///< (rad)
	float _altitude{NAN};   ///< above ellipsoid (m)
};
