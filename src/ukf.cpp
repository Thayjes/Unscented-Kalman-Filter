#include "ukf.h"
#include "measurement_package.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
    P_ = MatrixXd::Identity(5, 5);
    cout << P_ << endl;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    // Dimension of state x
    n_x_ = 5;
    // Augmented dimension of state with linear and yaw acceleration noise
    n_aug_ = 7;
    // Spreading parameter for sigma points
    lambda_ = 3 - n_aug_;
    
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    weights_ = VectorXd(2 * n_aug_ + 1);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        // velocity, yaw and yaw rate initial value
        x_(2) = 1; //
        x_(3) = 0; // Assuming straight heading
        x_(4) = 0; // this corresponds to  6 deg/sec
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float rho, phi, rho_dot;
            rho = measurement_pack.raw_measurements_[0];
            phi = measurement_pack.raw_measurements_[1];
            rho_dot = measurement_pack.raw_measurements_[2];
            x_(0) = rho*cos(phi);
            x_(1) = rho*sin(phi);
            float vx = rho_dot * cos(phi);
            float vy = rho_dot * sin(phi);
            float v  = sqrt(vx * vx + vy * vy);
            x_(2) = v;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            x_(0) = measurement_pack.raw_measurements_[0];
            x_(1) = measurement_pack.raw_measurements_[1];
        }
        time_us_ = measurement_pack.timestamp_;
        
        // set weights (this can be done here as they are same throughout)
        double weight_0 = lambda_/(lambda_+n_aug_);
        weights_(0) = weight_0;
        double weight = 0.5/(n_aug_+lambda_);
        for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
            weights_(i) = weight;
        }
        

        // done initializing, no need to predict or update
        is_initialized_ = true;
        cout << "Init " << endl;
        cout << x_ << endl;
        return;
    }
    
    
    
    ///* PREDICT THE SIGMA POINTS

    // Find the change in time between measurements
    double delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
    time_us_ = measurement_pack.timestamp_;
    cout << "Delta t = " << delta_t << endl;
    // Generate and Predict the sigma points
    Prediction(delta_t);
    cout << "predict:" << endl;
    cout << "x_" << endl << x_ << endl;

    
    // Update the state and covariance matrix based on the type of measurement
    if(measurement_pack.sensor_type_ == MeasurementPackage::LASER & use_laser_)
    {
        UpdateLidar(measurement_pack);
    }
    
    if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR & use_radar_)
    {
        UpdateRadar(measurement_pack);
    }
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    ///* GENERATE SIGMA POINTS
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    // Create augmented state vector
    Eigen::VectorXd x_aug(7);
    cout << "x_ before generation = " << endl << x_ << endl;
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    cout << "p aug = " << endl;
    cout << P_aug << endl;
    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    ///* END OF GENERATE SIGMA POINTS
    
    ///* PREDICT SIGMA POINTS
    //predict sigma points
    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
    ///* END OF PREDICT SIGMA POINTS
    
    ///* PREDICT MEAN AND COVARIANCE
    Xsig_pred_ = Xsig_pred;
    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);
    
    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);
    
    
    
    //predicted state mean
    x.fill(0.0);
    //predicted state covariance matrix
    P.fill(0.0);

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x = x+ weights_(i) * Xsig_pred_.col(i);
    }
    //predicted state covariance matrix
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P = P + weights_(i) * x_diff * x_diff.transpose() ;
    }
    // Update x_ and P_
    x_ = x;
    P_ = P;
    cout << "Predicted mean = " << endl << x_ << endl;
    cout << "Predicted Covariance = " << endl <<  P_ << endl;
    
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    //set measurement dimension, radar can measure px and py 
    int n_z = 2;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        //double v  = Xsig_pred_(2,i);
        //double yaw = Xsig_pred_(3,i);
        
        //double v1 = cos(yaw)*v;
        //double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = p_x;                       //x position
        Zsig(1,i) = p_y;                                //y position
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<  std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
    S = S + R;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    // Now calculate the RADAR NIS to see how good our uncertainty of the predicted measurement is
    // note that this is a combination of state and measurement uncertainty (Q and R)
    MatrixXd S_inv = S.inverse();
    nis_lidar_ = z_diff.transpose()*S_inv*z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        double rho = sqrt(p_x*p_x + p_y*p_y);
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        if(fabs(p_x) > 0.001)
            Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        if(fabs(rho) > 0.001)
            Zsig(2,i) = (p_x*v1 + p_y*v2 ) / rho;   //r_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    // Now calculate the RADAR NIS to see how good our uncertainty of the predicted measurement is
    // note that this is a combination of state and measurement uncertainty (Q and R)
    MatrixXd S_inv = S.inverse();
    nis_radar_ = z_diff.transpose()*S_inv*z_diff;
    
}
