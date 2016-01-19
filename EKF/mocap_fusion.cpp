/****************************************************************************
 *
 *   Copyright (c) 2015 Estimation and Control Library (ECL). All rights reserved.
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
 * 3. Neither the name ECL nor the names of its contributors may be
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

/**
 * @file mocap_fusion.cpp
 * Function for fusing motion capture (mocap) measurements/
 *
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Søren Andersen <s123369@student.dtu.dk>
 *
 */

#include "ekf.h"

void Ekf::fuseMocap()
{
	float innovations[6] = {};	// position in NED (x,y,z)->(E,N,D)
	float R[6] = {};			// sigma^2
	float Kfusion[24] = {};		// Array for storing calculated kalman gain

	// Calculate innovations
	innovations[0] = _state.pos(0) - _mocap_sample_delayed.position(0);
	innovations[1] = _state.pos(1) - _mocap_sample_delayed.position(1);
	innovations[2] = _state.pos(2) - _mocap_sample_delayed.position(2);
	R[0] = _params.mocap_x_noise;
	R[1] = _params.mocap_y_noise;
	R[2] = _params.mocap_z_noise;

	// Calculate quaternion error
	Quaternion quat_inv = _state.quat_nominal.inversed();
	Quaternion q_error =  _mocap_sample_delayed.attitude * quat_inv;
	q_error.normalize();
	Vector3f delta_ang_error;

	float scalar;

	if (q_error(0) >= 0.0f) {
		scalar = -2.0f;
		
	} else {
		scalar = 2.0f;
	}

	innovations[3] = scalar * q_error(1);
	innovations[4] = scalar * q_error(2);
	innovations[5] = scalar * q_error(3);
	R[3] = _params.mocap_r_noise;
	R[4] = _params.mocap_p_noise;
	R[5] = _params.mocap_h_noise;

	// Do fusion for each axis for position
	for (unsigned obs_index = 0; obs_index < 6; obs_index++) {

		unsigned state_index = obs_index<3?obs_index + 6:obs_index - 3;	// we start with x and this is the 7. state, then do r_error (state 0)

		// compute the innovation variance SK = HPH + R
		float S = P[state_index][state_index] + R[obs_index];
		S = 1.0f / S;

		// calculate kalman gain K = PHS
		for (int row = 0; row < 24; row++) {
			Kfusion[row] = P[row][state_index] * S;
		}

		// by definition the angle error state is zero at the fusion time
		_state.ang_error.setZero();

		// fuse the observation
		fuse(Kfusion, innovations[obs_index]);

		// correct the nominal quaternion
		Quaternion dq;
		dq.from_axis_angle(_state.ang_error);
		_state.quat_nominal = dq * _state.quat_nominal;
		_state.quat_nominal.normalize();

		// update covarinace matrix via Pnew = (I - KH)P
		float KHP[_k_num_states][_k_num_states] = {};

		for (unsigned row = 0; row < _k_num_states; row++) {
			for (unsigned column = 0; column < _k_num_states; column++) {
				KHP[row][column] = Kfusion[row] * P[state_index][column];
			}
		}

		for (unsigned row = 0; row < _k_num_states; row++) {
			for (unsigned column = 0; column < _k_num_states; column++) {
				P[row][column] = P[row][column] - KHP[row][column];
			}
		}

		makeSymmetrical();
		limitCov();
	}

}

// Ekf::fuse(float *K, float innovation) is reused from vel_pos_fusion
