#ifndef TRAJ_UTILS_HPP
#define TRAJ_UTILS_HPP

#include "root_finder.hpp"
#include <vector>
#include <list>
#include <Eigen/Eigen>

// Polynomial order and trajectory dimension are fixed here
constexpr int TrajOrder = 5;
constexpr int TrajDim = 3;

// Type for piece boundary condition and coefficient matrix
typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> BoundaryCond;
typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> CoefficientMat;
typedef Eigen::Matrix<double, TrajDim, TrajOrder> VelCoefficientMat;
typedef Eigen::Matrix<double, TrajDim, TrajOrder - 1> AccCoefficientMat;
typedef Eigen::Matrix<double, TrajDim, TrajOrder - 2> JerkCoefficientMat;

typedef Eigen::Matrix<double, 6, 1> StatePV;
typedef Eigen::Matrix<double, 9, 1> StatePVA;
typedef Eigen::Matrix<double, 10, 1> StatePVAM;
typedef Eigen::Matrix<double, TrajDim, 1> ControlJrk;
typedef Eigen::Matrix<double, TrajDim, 1> ControlAcc;

// A single piece of a trajectory, which is indeed a polynomial
class Piece
{
private:
    // Piece(t) = c5*t^5 + c4*t^4 + ... + c1*t + c0
    // The natural coefficient matrix = [c5,c4,c3,c2,c1,c0]
    double duration;
    // Any time in [0, T] is normalized into [0.0, 1.0]
    // Therefore, nCoeffMat = [c5*T^5,c4*T^4,c3*T^3,c2*T^2,c1*T,c0*1]
    // is used for better numerical stability
    CoefficientMat nCoeffMat;

public:
    Piece() = default;

    // Constructor from duration and coefficient
    Piece(double dur, const CoefficientMat &coeffs) : duration(dur)
    {
        double t = 1.0;
        for (int i = TrajOrder; i >= 0; i--)
        {
            nCoeffMat.col(i) = coeffs.col(i) * t;
            t *= dur;
        }
    }

    // Constructor from boundary condition and duration
    Piece(const BoundaryCond &boundCond, double dur) : duration(dur)
    {
        // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
        double t1 = dur;
        double t2 = t1 * t1;

        // Inverse mapping is computed without explicit matrix inverse
        // It maps boundary condition to normalized coefficient matrix
        nCoeffMat.col(0) = 0.5 * (boundCond.col(5) - boundCond.col(2)) * t2 -
                           3.0 * (boundCond.col(1) + boundCond.col(4)) * t1 +
                           6.0 * (boundCond.col(3) - boundCond.col(0));
        nCoeffMat.col(1) = (-boundCond.col(5) + 1.5 * boundCond.col(2)) * t2 +
                           (8.0 * boundCond.col(1) + 7.0 * boundCond.col(4)) * t1 +
                           15.0 * (-boundCond.col(3) + boundCond.col(0));
        nCoeffMat.col(2) = (0.5 * boundCond.col(5) - 1.5 * boundCond.col(2)) * t2 -
                           (6.0 * boundCond.col(1) + 4.0 * boundCond.col(4)) * t1 +
                           10.0 * (boundCond.col(3) - boundCond.col(0));
        nCoeffMat.col(3) = 0.5 * boundCond.col(2) * t2;
        nCoeffMat.col(4) = boundCond.col(1) * t1;
        nCoeffMat.col(5) = boundCond.col(0);
    }

    inline int getDim() const
    {
        return TrajDim;
    }

    inline int getOrder() const
    {
        return TrajOrder;
    }

    inline double getDuration() const
    {
        return duration;
    }

    // Get the position at time t in this piece
    inline Eigen::Vector3d getPos(double t) const
    {
        // Normalize the time
        t /= duration;
        Eigen::Vector3d pos(0.0, 0.0, 0.0);
        double tn = 1.0;
        for (int i = TrajOrder; i >= 0; i--)
        {
            pos += tn * nCoeffMat.col(i);
            tn *= t;
        }
        // The pos is not affected by normalization
        return pos;
    }

    // Get the velocity at time t in this piece
    inline Eigen::Vector3d getVel(double t) const
    {
        // Normalize the time
        t /= duration;
        Eigen::Vector3d vel(0.0, 0.0, 0.0);
        double tn = 1.0;
        int n = 1;
        for (int i = TrajOrder - 1; i >= 0; i--)
        {
            vel += n * tn * nCoeffMat.col(i);
            tn *= t;
            n++;
        }
        // Recover the actual vel
        vel /= duration;
        return vel;
    }

    // Get the acceleration at time t in this piece
    inline Eigen::Vector3d getAcc(double t) const
    {
        // Normalize the time
        t /= duration;
        Eigen::Vector3d acc(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        for (int i = TrajOrder - 2; i >= 0; i--)
        {
            acc += m * n * tn * nCoeffMat.col(i);
            tn *= t;
            m++;
            n++;
        }
        // Recover the actual acc
        acc /= duration * duration;
        return acc;
    }

    // Get the jerk at time t in this piece
    inline Eigen::Vector3d getJerk(double t) const
    {
        // Normalize the time
        t /= duration;
        Eigen::Vector3d jerk(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        int k = 3;
        for (int i = TrajOrder - 3; i >= 0; i--)
        {
            jerk += k * m * n * tn * nCoeffMat.col(i);
            tn *= t;
            k++;
            m++;
            n++;
        }
        // Recover the actual acc
        jerk /= duration * duration * duration;
        return jerk;
    }

    // Get the boundary condition of this piece
    inline BoundaryCond getBoundCond() const
    {
        BoundaryCond boundCond;
        boundCond << getPos(0.0), getVel(0.0), getAcc(0.0),
            getPos(duration), getVel(duration), getAcc(duration);
        return boundCond;
    }

    // Get the coefficient matrix of the piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline CoefficientMat getCoeffMat(bool normalized = false) const
    {
        CoefficientMat posCoeffsMat;
        double t = 1;
        for (int i = TrajOrder; i >= 0; i--)
        {
            posCoeffsMat.col(i) = nCoeffMat.col(i) / t;
            t *= normalized ? 1.0 : duration;
        }
        return posCoeffsMat;
    }

    // Get the polynomial coefficients of velocity of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline VelCoefficientMat getVelCoeffMat(bool normalized = false) const
    {
        VelCoefficientMat velCoeffMat;
        int n = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration;
        for (int i = TrajOrder - 1; i >= 0; i--)
        {
            velCoeffMat.col(i) = n * nCoeffMat.col(i) / t;
            n++;
            t *= normalized ? 1.0 : duration;
        }
        return velCoeffMat;
    }

    // Get the polynomial coefficients of acceleration of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline AccCoefficientMat getAccCoeffMat(bool normalized = false) const
    {
        AccCoefficientMat accCoeffMat;
        int n = 2;
        int m = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration * duration;
        for (int i = TrajOrder - 2; i >= 0; i--)
        {
            accCoeffMat.col(i) = n * m * nCoeffMat.col(i) / t;
            n++;
            m++;
            t *= normalized ? 1.0 : duration;
        }
        return accCoeffMat;
    }

    // Get the polynomial coefficients of jerk of this piece
    // Default arg chooses the natural coefficients
    // If normalized version is needed, set the arg true
    inline JerkCoefficientMat getJerkCoeffMat(bool normalized = false) const
    {
        JerkCoefficientMat jerkCoeffMat;
        int n = 3;
        int m = 2;
        int l = 1;
        double t = 1.0;
        t *= normalized ? 1.0 : duration * duration * duration;
        for (int i = TrajOrder - 3; i >= 0; i--)
        {
            jerkCoeffMat.col(i) = n * m * l * nCoeffMat.col(i) / t;
            n++;
            m++;
            l++;
            t *= normalized ? 1.0 : duration;
        }
        return jerkCoeffMat;
    }

    // Get the max velocity rate of the piece
    inline double getMaxVelRate() const
    {
        // Compute normalized squared vel norm polynomial coefficient matrix
        Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++)
        {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
        {
            return 0.0;
        }
        else
        {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
            {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
            {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);

            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxVelRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++)
            {
                if (0.0 <= *it && 1.0 >= *it)
                {
                    // Recover the actual time then get the vel squared norm
                    tempNormSqr = getVel((*it) * duration).squaredNorm();
                    maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
                }
            }
            return sqrt(maxVelRateSqr);
        }
    }

    // Get the max acceleration rate of the piece
    inline double getMaxAccRate() const
    {
        // Compute normalized squared acc norm polynomial coefficient matrix
        Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                RootFinder::polySqr(nAccCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++)
        {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
        {
            return 0.0;
        }
        else
        {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
            {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
            {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxAccRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++)
            {
                if (0.0 <= *it && 1.0 >= *it)
                {
                    // Recover the actual time then get the acc squared norm
                    tempNormSqr = getAcc((*it) * duration).squaredNorm();
                    maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
                }
            }
            return sqrt(maxAccRateSqr);
        }
    }

    inline double getMaxJerkRate() const
    {
        // Compute normalized squared jerk norm polynomial coefficient matrix
        Eigen::MatrixXd nJerkCoeffMat = getJerkCoeffMat(true);
        Eigen::VectorXd coeff = RootFinder::polySqr(nJerkCoeffMat.row(0)) +
                                RootFinder::polySqr(nJerkCoeffMat.row(1)) +
                                RootFinder::polySqr(nJerkCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++)
        {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
        {
            return 0.0;
        }
        else
        {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
            {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
            {
                r = 0.5 * (r + 1.0);
            }
            // Find all stationaries
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            // Check boundary points and stationaries within duration
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxJerkRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++)
            {
                if (0.0 <= *it && 1.0 >= *it)
                {
                    // Recover the actual time then get the acc squared norm
                    tempNormSqr = getJerk((*it) * duration).squaredNorm();
                    maxJerkRateSqr = maxJerkRateSqr < tempNormSqr ? tempNormSqr : maxJerkRateSqr;
                }
            }
            return sqrt(maxJerkRateSqr);
        }
    }

    // Check whether velocity rate of the piece is always less than maxVelRate
    inline bool checkMaxVelRate(double maxVelRate) const
    {
        double sqrMaxVelRate = maxVelRate * maxVelRate;
        if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
            getVel(duration).squaredNorm() >= sqrMaxVelRate)
        {
            return false;
        }
        else
        {
            Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(2));
            // Convert the actual squared maxVelRate to a normalized one
            double t2 = duration * duration;
            coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    // Check whether accleration rate of the piece is always less than maxAccRate
    inline bool checkMaxAccRate(double maxAccRate) const
    {
        double sqrMaxAccRate = maxAccRate * maxAccRate;
        if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
            getAcc(duration).squaredNorm() >= sqrMaxAccRate)
        {
            return false;
        }
        else
        {
            Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(2));
            // Convert the actual squared maxAccRate to a normalized one
            double t2 = duration * duration;
            double t4 = t2 * t2;
            coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    // Check whether jerk rate of the piece is always less than maxJerkRate
    inline bool checkMaxJerkRate(double maxJerkRate) const
    {
        double sqrMaxJerkRate = maxJerkRate * maxJerkRate;
        if (getJerk(0.0).squaredNorm() >= sqrMaxJerkRate ||
            getJerk(duration).squaredNorm() >= sqrMaxJerkRate)
        {
            return false;
        }
        else
        {
            Eigen::MatrixXd nJerkCoeffMat = getJerkCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nJerkCoeffMat.row(0)) +
                                    RootFinder::polySqr(nJerkCoeffMat.row(1)) +
                                    RootFinder::polySqr(nJerkCoeffMat.row(2));
            // Convert the actual squared maxJerkRate to a normalized one
            double t2 = duration * duration;
            double t4 = t2 * t2;
            double t6 = t4 * t2;
            coeff.tail<1>()(0) -= sqrMaxJerkRate * t6;
            // Directly check the root existence in the normalized interval
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    //Scale the Piece(t) to Piece(k*t)
    inline void scaleTime(double k)
    {
        duration /= k;
        return;
    }

    inline void sampleOneSeg(std::vector< StatePVA >* vis_x) const 
    {
        double dt = 0.005;
        for (double t = 0.0; t < duration; t += dt) 
        {
            Eigen::Vector3d pos, vel, acc;
            pos = getPos(t);
            vel = getVel(t);
            acc = getAcc(t);
            StatePVA x;
            x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2);
            vis_x->push_back(x);
        }
    }

    inline void sampleOneSeg(std::vector< StatePVAM >* vis_x, double sample_t, double ground_judge) const 
    {
        double dt = sample_t;
        for (double t = 0.0; t < duration; t += dt) 
        {
            Eigen::Vector3d pos, vel, acc;
            pos = getPos(t);
            vel = getVel(t);
            acc = getAcc(t);
            int motion_state = 0;
            if(pos[2] >= ground_judge) motion_state = 1;
            StatePVAM x;
            x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2), motion_state;
            vis_x->push_back(x);
        }
    }

    // for 5th degree polynomial
    inline void cutPiece(const Piece &orig_piece, double ts, CoefficientMat &new_coeff) const
    {
        CoefficientMat ori_coeff = orig_piece.getCoeffMat();
        double ts2 = ts * ts;
        double ts3 = ts2 * ts;
        double ts4 = ts3 * ts;
        double ts5 = ts4 * ts;
        for (int dim = 0; dim < 3; ++dim)
        {
            new_coeff(dim, 0) = ori_coeff(dim, 0);  //c5*t^5
            new_coeff(dim, 1) = ori_coeff(dim, 1) + 5*ori_coeff(dim, 0)*ts;  //c4*4^4
            new_coeff(dim, 2) = ori_coeff(dim, 2) + 4*ori_coeff(dim, 1)*ts + 10*ori_coeff(dim, 0)*ts2;
            new_coeff(dim, 3) = ori_coeff(dim, 3) + 3*ori_coeff(dim, 2)*ts + 6*ori_coeff(dim, 1)*ts2 + 10*ori_coeff(dim, 0)*ts3;
            new_coeff(dim, 4) = ori_coeff(dim, 4) + 2*ori_coeff(dim, 3)*ts + 3*ori_coeff(dim, 2)*ts2 + 4*ori_coeff(dim, 1)*ts3 + 5*ori_coeff(dim, 0)*ts4;
            new_coeff(dim, 5) = ori_coeff(dim, 5) + ori_coeff(dim, 4)*ts + ori_coeff(dim, 3)*ts2 + ori_coeff(dim, 2)*ts3 + ori_coeff(dim, 1)*ts4 + ori_coeff(dim, 0)*ts5;
        }
    }

    // for 5th degree polynomial， cost = integral(j^T rho j) + tau
    inline double calCost(const double &rho) const
    {
        double tau2 = duration * duration;
        double tau3 = tau2 * duration;
        double tau4 = tau3 * duration;
        double tau5 = tau4 * duration;

        CoefficientMat coeff = getCoeffMat();
        Eigen::Matrix<double, 6, 6> B = Eigen::Matrix<double, 6, 6>::Zero(6, 6);
        B(0, 0) = 720 * tau5;
        B(1, 1) = 192 * tau3;
        B(2, 2) = 36 * duration;
        B(0, 1) = B(1, 0) = 360 * tau4;
        B(0, 2) = B(2, 0) = 120 * tau3;
        B(1, 2) = B(2, 1) = 72 * tau2;
        double cost(0.0);
        for (int i=0; i<3; i++)
        {
            cost += coeff.row(i) * B * coeff.row(i).transpose();
        }
        cost *= rho;
        cost += duration;

        return cost;
    }

    inline double project_pt(const Eigen::Vector3d &pt,
                           double &tt, Eigen::Vector3d &pro_pt) {
        // 2*(p-p0)^T * \dot{p} = 0
        auto l_coeff = getCoeffMat();
        l_coeff.col(5) = l_coeff.col(5) - pt;
        auto r_coeff = getVelCoeffMat();
        Eigen::VectorXd eq = Eigen::VectorXd::Zero(2 * 5);
        for (int j = 0; j < l_coeff.rows(); ++j) {
            eq = eq + RootFinder::polyConv(l_coeff.row(j), r_coeff.row(j));
        }
        double l = -0.0625;
        double r = duration + 0.0625;
        while (fabs(RootFinder::polyVal(eq, l)) < DBL_EPSILON) {
            l = 0.5 * l;
        }
        while (fabs(RootFinder::polyVal(eq, r)) < DBL_EPSILON) {
            r = 0.5 * (duration + r);
        }
        std::set<double> roots =
            RootFinder::solvePolynomial(eq, l, r, 1e-6);
        // std::cout << "# roots: " << roots.size() << std::endl;
        double min_dist = -1;
        for (const auto &root : roots) {
            // std::cout << "root: " << root << std::endl;
            if (root < 0 || root > duration) {
                continue;
            }
            if (getVel(root).norm() < 1e-6) { // velocity == 0, ignore it
                continue;
            }
            // std::cout << "find min!" << std::endl;
            Eigen::Vector3d p = getPos(root);
            // std::cout << "p: " << p.transpose() << std::endl;
            double distance = (p - pt).norm();
            if (distance < min_dist || min_dist < 0) {
                min_dist = distance;
                tt = root;
                pro_pt = p;
            }
        }
        return min_dist;
    }
    inline bool intersection_plane(const Eigen::Vector3d p, 
                                   const Eigen::Vector3d v,
                                   double &tt, Eigen::Vector3d &pt) const {
        // (pt - p)^T * v = 0
        auto coeff = getCoeffMat();
        coeff.col(5) = coeff.col(5) - p;
        Eigen::VectorXd eq = coeff.transpose() * v;
        double l = -0.0625;
        double r = duration + 0.0625;
        while (fabs(RootFinder::polyVal(eq, l)) < DBL_EPSILON) {
            l = 0.5 * l;
        }
        while (fabs(RootFinder::polyVal(eq, r)) < DBL_EPSILON) {
            r = 0.5 * (duration + r);
        }
        std::set<double> roots =
            RootFinder::solvePolynomial(eq, l, r, 1e-6);
        for (const auto &root : roots) {
            tt = root;
            pt = getPos(root);
            return true;
        }
        return false;
    }

};

// A whole trajectory which contains multiple pieces
class Trajectory
{
private:
    typedef std::vector<Piece> Pieces;
    Pieces pieces;

public:
    Trajectory() = default;

    // Constructor from durations and coefficient matrices
    Trajectory(const std::vector<double> &durs,
               const std::vector<CoefficientMat> &coeffMats)
    {
        int N = std::min(durs.size(), coeffMats.size());
        pieces.reserve(N);
        for (int i = 0; i < N; i++)
        {
            pieces.emplace_back(durs[i], coeffMats[i]);
        }
    }

    inline int getPieceNum() const
    {
        return pieces.size();
    }

    // Get durations vector of all pieces
    inline std::vector<double> getDurations() const
    {
        std::vector<double> durations;
        durations.reserve(getPieceNum());
        for (int i = 0; i < getPieceNum(); i++)
        {
            durations.push_back(pieces[i].getDuration());
        }
        return durations;
    }

    // Get total duration of the trajectory
    inline double getTotalDuration() const
    {
        double totalDuration = 0.0;
        for (int i = 0; i < getPieceNum(); i++)
        {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }

    // Reload the operator[] to reach the i-th piece
    inline const Piece &operator[](int i) const
    {
        return pieces[i];
    }

    inline Piece &operator[](int i)
    {
        return pieces[i];
    }

    inline void clear(void)
    {
        pieces.clear();
    }

    inline Pieces::const_iterator begin() const
    {
        return pieces.begin();
    }

    inline Pieces::const_iterator end() const
    {
        return pieces.end();
    }

    inline void reserve(const int &n)
    {
        pieces.reserve(n);
        return;
    }

    // Put another piece at the tail of this trajectory
    inline void emplace_back(const Piece &piece)
    {
        pieces.emplace_back(piece);
        return;
    }

    // Two corresponding constructors of Piece both are supported here
    template <typename ArgTypeL, typename ArgTypeR>
    inline void emplace_back(const ArgTypeL &argL, const ArgTypeR &argR)
    {
        pieces.emplace_back(argL, argR);
        return;
    }

    // Append another Trajectory at the tail of this trajectory
    inline void append(const Trajectory &traj)
    {
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }

    // Find the piece at which the time t is located
    // The index is returned and the offset in t is removed
    inline int locatePieceIdx(double &t) const
    {
        int idx;
        double dur;
        for (idx = 0;
             idx < getPieceNum() &&
             t > (dur = pieces[idx].getDuration());
             idx++)
        {
            t -= dur;
        }
        if (idx == getPieceNum())
        {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    // Get the position at time t of the trajectory
    inline Eigen::Vector3d getPos(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    // Get the velocity at time t of the trajectory
    inline Eigen::Vector3d getVel(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getVel(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getAcc(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getAcc(t);
    }

    // Get the position at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncPos(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getPos(0.0);
        }
        else
        {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the velocity at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncVel(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getVel(0.0);
        }
        else
        {
            return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the acceleration at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncAcc(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getAcc(0.0);
        }
        else
        {
            return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the max velocity rate of the trajectory
    inline double getMaxVelRate() const
    {
        double maxVelRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++)
        {
            tempNorm = pieces[i].getMaxVelRate();
            maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
        }
        return maxVelRate;
    }

    // Get the max acceleration rate of the trajectory
    inline double getMaxAccRate() const
    {
        double maxAccRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++)
        {
            tempNorm = pieces[i].getMaxAccRate();
            maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
        }
        return maxAccRate;
    }

    // Get the max jerk rate of the trajectory
    inline double getMaxJerkRate() const
    {
        double maxJerkRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++)
        {
            tempNorm = pieces[i].getMaxJerkRate();
            maxJerkRate = maxJerkRate < tempNorm ? tempNorm : maxJerkRate;
        }
        return maxJerkRate;
    }

    // Check whether the velocity rate of this trajectory exceeds the threshold
    inline bool checkMaxVelRate(double maxVelRate) const
    {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++)
        {
            feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
        }
        return feasible;
    }

    // Check whether the acceleration rate of this trajectory exceeds the threshold
    inline bool checkMaxAccRate(double maxAccRate) const
    {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++)
        {
            feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
        }
        return feasible;
    }

    // Check whether the jerk rate of this trajectory exceeds the threshold
    inline bool checkMaxJerkRate(double maxJerkRate) const
    {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++)
        {
            feasible = feasible && pieces[i].checkMaxJerkRate(maxJerkRate);
        }
        return feasible;
    }

    // Scale the Trajectory(t) to Trajectory(k*t)
    inline void scaleTime(double k)
    {
        for (int i = 0; i < getPieceNum(); i++)
        {
            pieces[i].scaleTime(k);
        }
    }

    inline void sampleWholeTrajectory(std::vector< StatePVA >* vis_x) const 
    {
        int n = getPieceNum();
        for (int i = 0; i < n; ++i)
        {
            pieces[i].sampleOneSeg(vis_x);
        }
    }

    inline double calCost(const double &rho, double* seg_cost) const
    {
        double cost(0.0);
        for (int i = 0; i < getPieceNum(); i++)
        {
            seg_cost[i] = pieces[i].calCost(rho);
            cost += seg_cost[i];
        }
        return cost;
    }

    inline void sampleWholeTrajectoryForOptimization(std::vector< StatePVAM >* vis_x, double sample_t, double ground_judge) const 
    {
        int n = getPieceNum();
        for (int i = 0; i < n; ++i)
        {
            pieces[i].sampleOneSeg(vis_x, sample_t, ground_judge);
        }
    }

    inline void getWpts(std::vector< StatePVA >* wpts)
    {
        Eigen::Vector3d pos, vel, acc;
        StatePVA x;
        pos = pieces[0].getPos(0);
        vel = pieces[0].getVel(0);
        acc = pieces[0].getAcc(0);
        x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2);
        wpts->push_back(x);
        
        int n = getPieceNum();
        for (int i = 0; i < n; ++i)
        {
            double t = pieces[i].getDuration();
            pos = pieces[i].getPos(t);
            vel = pieces[i].getVel(t);
            acc = pieces[i].getAcc(t);
            x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2);
            wpts->push_back(x);
        }
    }

    inline const Piece& getPiece(int i) const {
        return pieces[i];
    }
    inline double project_pt(const Eigen::Vector3d &pt,
                           int &ii, double &tt, Eigen::Vector3d &pro_pt) {
        double dist = -1;
        for (int i=0; i<getPieceNum(); ++i) {
            auto piece = pieces[i];
            dist = piece.project_pt(pt, tt, pro_pt);
            if (dist > 0) {
                ii = i;
                break;
            }
        }
        // if (dist < 0) {
        //     std::cout << "\033[32m" << "cannot project pt to traj" << "\033[0m" << std::endl;
        //     // std::cout << "pt: " << pt.transpose() << std::endl;
        //     // assert(false);
        // }
        return dist;
    }
    inline bool intersection_plane(const Eigen::Vector3d p, 
                                   const Eigen::Vector3d v,
                                   int &ii, double &tt, Eigen::Vector3d &pt) {
        for (int i=0; i<getPieceNum(); ++i) {
            const auto& piece = pieces[i];
            if ( piece.intersection_plane(p,v,tt,pt) ) {
                ii = i;
                return true;
            }
        }
        return false;
    }

    inline double evaluateTrajJerk() const
    {
        double objective = 0.0;
        int M = getPieceNum();
        CoefficientMat cMat;
        double t1, t2, t3, t4, t5;
        for (int i = 0; i < M; i++)
        {
            cMat = operator[](i).getCoeffMat();
            t1 = operator[](i).getDuration();
            t2 = t1 * t1;
            t3 = t2 * t1;
            t4 = t2 * t2;
            t5 = t2 * t3;
            objective += 36.0 * cMat.col(2).squaredNorm() * t1 +
                        144.0 * cMat.col(1).dot(cMat.col(2)) * t2 +
                        192.0 * cMat.col(1).squaredNorm() * t3 +
                        240.0 * cMat.col(0).dot(cMat.col(2)) * t3 +
                        720.0 * cMat.col(0).dot(cMat.col(1)) * t4 +
                        720.0 * cMat.col(0).squaredNorm() * t5;
        }
        return objective;
    }
};

// The banded system class is used for solving
// banded linear system Ax=b efficiently.
// A is an N*N band matrix with lower band width lowerBw
// and upper band width upperBw.
// Banded LU factorization has O(N) time complexity.
class BandedSystem
{
public:
    // The size of A, as well as the lower/upper
    // banded width p/q are needed
    inline void create(const int &n, const int &p, const int &q)
    {
        // In case of re-creating before destroying
        destroy();
        N = n;
        lowerBw = p;
        upperBw = q;
        int actualSize = N * (lowerBw + upperBw + 1);
        ptrData = new double[actualSize];
        std::fill_n(ptrData, actualSize, 0.0);
        return;
    }

    inline void destroy()
    {
        if (ptrData != nullptr)
        {
            delete[] ptrData;
            ptrData = nullptr;
        }
        return;
    }

private:
    int N;
    int lowerBw;
    int upperBw;
    // Compulsory nullptr initialization here
    double *ptrData = nullptr;

public:
    // Reset the matrix to zero
    inline void reset(void)
    {
        std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
        return;
    }

    // The band matrix is stored as suggested in "Matrix Computation"
    inline const double &operator()(const int &i, const int &j) const
    {
        return ptrData[(i - j + upperBw) * N + j];
    }

    inline double &operator()(const int &i, const int &j)
    {
        return ptrData[(i - j + upperBw) * N + j];
    }

    // This function conducts banded LU factorization in place
    // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
    inline void factorizeLU()
    {
        int iM, jM;
        double cVl;
        for (int k = 0; k <= N - 2; k++)
        {
            iM = std::min(k + lowerBw, N - 1);
            cVl = operator()(k, k);
            for (int i = k + 1; i <= iM; i++)
            {
                if (operator()(i, k) != 0.0)
                {
                    operator()(i, k) /= cVl;
                }
            }
            jM = std::min(k + upperBw, N - 1);
            for (int j = k + 1; j <= jM; j++)
            {
                cVl = operator()(k, j);
                if (cVl != 0.0)
                {
                    for (int i = k + 1; i <= iM; i++)
                    {
                        if (operator()(i, k) != 0.0)
                        {
                            operator()(i, j) -= operator()(i, k) * cVl;
                        }
                    }
                }
            }
        }
        return;
    }

    // This function solves Ax=b, then stores x in b
    // The input b is required to be N*m, i.e.,
    // m vectors to be solved.
    inline void solve(Eigen::MatrixXd &b) const
    {
        int iM;
        for (int j = 0; j <= N - 1; j++)
        {
            iM = std::min(j + lowerBw, N - 1);
            for (int i = j + 1; i <= iM; i++)
            {
                if (operator()(i, j) != 0.0)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; j--)
        {
            b.row(j) /= operator()(j, j);
            iM = std::max(0, j - upperBw);
            for (int i = iM; i <= j - 1; i++)
            {
                if (operator()(i, j) != 0.0)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        return;
    }

    // This function solves ATx=b, then stores x in b
    // The input b is required to be N*m, i.e.,
    // m vectors to be solved.
    inline void solveAdj(Eigen::MatrixXd &b) const
    {
        int iM;
        for (int j = 0; j <= N - 1; j++)
        {
            b.row(j) /= operator()(j, j);
            iM = std::min(j + upperBw, N - 1);
            for (int i = j + 1; i <= iM; i++)
            {
                if (operator()(j, i) != 0.0)
                {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; j--)
        {
            iM = std::max(0, j - lowerBw);
            for (int i = iM; i <= j - 1; i++)
            {
                if (operator()(j, i) != 0.0)
                {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        return;
    }
};

class MinJerkOpt
{
public:
    MinJerkOpt() = default;
    ~MinJerkOpt() { A.destroy(); }

private:
    int N;
    Eigen::Matrix3d headPVA;
    Eigen::Matrix3d tailPVA;
    Eigen::VectorXd T1;
    BandedSystem A;
    Eigen::MatrixXd b;

    // Temp variables
    Eigen::VectorXd T2;
    Eigen::VectorXd T3;
    Eigen::VectorXd T4;
    Eigen::VectorXd T5;
    Eigen::MatrixXd gdC;

private:
    template <typename EIGENVEC>
    inline void addGradJbyT(EIGENVEC &gdT) const
    {
        for (int i = 0; i < N; i++)
        {
            gdT(i) += 36.0 * b.row(6 * i + 3).squaredNorm() +
                      288.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T1(i) +
                      576.0 * b.row(6 * i + 4).squaredNorm() * T2(i) +
                      720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T2(i) +
                      2880.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T3(i) +
                      3600.0 * b.row(6 * i + 5).squaredNorm() * T4(i);
        }
        return;
    }

    template <typename EIGENMAT>
    inline void addGradJbyC(EIGENMAT &gdC) const
    {
        for (int i = 0; i < N; i++)
        {
            gdC.row(6 * i + 5) += 240.0 * b.row(6 * i + 3) * T3(i) +
                                  720.0 * b.row(6 * i + 4) * T4(i) +
                                  1440.0 * b.row(6 * i + 5) * T5(i);
            gdC.row(6 * i + 4) += 144.0 * b.row(6 * i + 3) * T2(i) +
                                  384.0 * b.row(6 * i + 4) * T3(i) +
                                  720.0 * b.row(6 * i + 5) * T4(i);
            gdC.row(6 * i + 3) += 72.0 * b.row(6 * i + 3) * T1(i) +
                                  144.0 * b.row(6 * i + 4) * T2(i) +
                                  240.0 * b.row(6 * i + 5) * T3(i);
        }
        return;
    }

    inline void solveAdjGradC(Eigen::MatrixXd &gdC) const
    {
        A.solveAdj(gdC);
        return;
    }

    template <typename EIGENVEC>
    inline void addPropCtoT(const Eigen::MatrixXd &adjGdC, EIGENVEC &gdT) const
    {
        Eigen::MatrixXd B1(6, 3), B2(3, 3);

        Eigen::RowVector3d negVel, negAcc, negJer, negSnp, negCrk;

        for (int i = 0; i < N - 1; i++)
        {
            negVel = -(b.row(i * 6 + 1) +
                       2.0 * T1(i) * b.row(i * 6 + 2) +
                       3.0 * T2(i) * b.row(i * 6 + 3) +
                       4.0 * T3(i) * b.row(i * 6 + 4) +
                       5.0 * T4(i) * b.row(i * 6 + 5));
            negAcc = -(2.0 * b.row(i * 6 + 2) +
                       6.0 * T1(i) * b.row(i * 6 + 3) +
                       12.0 * T2(i) * b.row(i * 6 + 4) +
                       20.0 * T3(i) * b.row(i * 6 + 5));
            negJer = -(6.0 * b.row(i * 6 + 3) +
                       24.0 * T1(i) * b.row(i * 6 + 4) +
                       60.0 * T2(i) * b.row(i * 6 + 5));
            negSnp = -(24.0 * b.row(i * 6 + 4) +
                       120.0 * T1(i) * b.row(i * 6 + 5));
            negCrk = -120.0 * b.row(i * 6 + 5);

            B1 << negSnp, negCrk, negVel, negVel, negAcc, negJer;

            gdT(i) += B1.cwiseProduct(adjGdC.block<6, 3>(6 * i + 3, 0)).sum();
        }

        negVel = -(b.row(6 * N - 5) +
                   2.0 * T1(N - 1) * b.row(6 * N - 4) +
                   3.0 * T2(N - 1) * b.row(6 * N - 3) +
                   4.0 * T3(N - 1) * b.row(6 * N - 2) +
                   5.0 * T4(N - 1) * b.row(6 * N - 1));
        negAcc = -(2.0 * b.row(6 * N - 4) +
                   6.0 * T1(N - 1) * b.row(6 * N - 3) +
                   12.0 * T2(N - 1) * b.row(6 * N - 2) +
                   20.0 * T3(N - 1) * b.row(6 * N - 1));
        negJer = -(6.0 * b.row(6 * N - 3) +
                   24.0 * T1(N - 1) * b.row(6 * N - 2) +
                   60.0 * T2(N - 1) * b.row(6 * N - 1));

        B2 << negVel, negAcc, negJer;

        gdT(N - 1) += B2.cwiseProduct(adjGdC.block<3, 3>(6 * N - 3, 0)).sum();

        return;
    }

    template <typename EIGENMAT>
    inline void addPropCtoP(const Eigen::MatrixXd &adjGdC, EIGENMAT &gdInP) const
    {
        for (int i = 0; i < N - 1; i++)
        {
            gdInP.col(i) += adjGdC.row(6 * i + 5).transpose();
        }
        return;
    }

    template <typename EIGENVEC>
    inline void addTimeIntPenalty(const Eigen::VectorXi cons,
                                  const Eigen::VectorXi &idxHs,
                                  const std::vector<Eigen::MatrixXd> &cfgHs,
                                  const double vmax,
                                  const double amax,
                                  const Eigen::Vector3d ci,
                                  double &cost,
                                  EIGENVEC &gdT,
                                  Eigen::MatrixXd &gdC) const
    {
        double pena = 0.0;
        const double vmaxSqr = vmax * vmax;
        const double amaxSqr = amax * amax;

        Eigen::Vector3d pos, vel, acc, jer;
        double step, alpha;
        double s1, s2, s3, s4, s5;
        Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3;
        Eigen::Vector3d outerNormal;
        int K;
        double violaPos, violaVel, violaAcc;
        double violaPosPenaD, violaVelPenaD, violaAccPenaD;
        double violaPosPena, violaVelPena, violaAccPena;
        Eigen::Matrix<double, 6, 3> gradViolaVc, gradViolaAc;
        double gradViolaVt, gradViolaAt;
        double omg;

        int innerLoop, idx;
        for (int i = 0; i < N; i++)
        {
            const auto &c = b.block<6, 3>(i * 6, 0);
            step = T1(i) / cons(i);
            s1 = 0.0;
            innerLoop = cons(i) + 1;
            for (int j = 0; j < innerLoop; j++)
            {
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                beta0 << 1.0, s1, s2, s3, s4, s5;
                beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
                beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
                beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
                alpha = 1.0 / cons(i) * j;
                pos = c.transpose() * beta0;
                vel = c.transpose() * beta1;
                acc = c.transpose() * beta2;
                jer = c.transpose() * beta3;
                violaVel = vel.squaredNorm() - vmaxSqr;
                violaAcc = acc.squaredNorm() - amaxSqr;

                omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;

                idx = idxHs(i);
                K = cfgHs[idx].cols();
                for (int k = 0; k < K; k++)
                {
                    outerNormal = cfgHs[idx].col(k).head<3>();
                    violaPos = outerNormal.dot(pos - cfgHs[idx].col(k).tail<3>());
                    if (violaPos > 0.0)
                    {
                        violaPosPenaD = violaPos * violaPos;
                        violaPosPena = violaPosPenaD * violaPos;
                        violaPosPenaD *= 3.0;
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(0) * violaPosPenaD * beta0 * outerNormal.transpose();
                        gdT(i) += omg * (ci(0) * violaPosPenaD * alpha * outerNormal.dot(vel) * step +
                                         ci(0) * violaPosPena / cons(i));
                        pena += omg * step * ci(0) * violaPosPena;
                    }
                }

                if (violaVel > 0.0)
                {
                    violaVelPenaD = violaVel * violaVel;
                    violaVelPena = violaVelPenaD * violaVel;
                    violaVelPenaD *= 3.0;
                    gradViolaVc = 2.0 * beta1 * vel.transpose();
                    gradViolaVt = 2.0 * alpha * vel.transpose() * acc;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * ci(1) * violaVelPenaD * gradViolaVc;
                    gdT(i) += omg * (ci(1) * violaVelPenaD * gradViolaVt * step +
                                     ci(1) * violaVelPena / cons(i));
                    pena += omg * step * ci(1) * violaVelPena;
                }

                if (violaAcc > 0.0)
                {
                    violaAccPenaD = violaAcc * violaAcc;
                    violaAccPena = violaAccPenaD * violaAcc;
                    violaAccPenaD *= 3.0;
                    gradViolaAc = 2.0 * beta2 * acc.transpose();
                    gradViolaAt = 2.0 * alpha * acc.transpose() * jer;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * ci(2) * violaAccPenaD * gradViolaAc;
                    gdT(i) += omg * (ci(2) * violaAccPenaD * gradViolaAt * step +
                                     ci(2) * violaAccPena / cons(i));
                    pena += omg * step * ci(2) * violaAccPena;
                }

                s1 += step;
            }
        }

        cost += pena;
        return;
    }

public:
    inline void reset(const Eigen::Matrix3d &headState,
                      const Eigen::Matrix3d &tailState,
                      const int &pieceNum)
    {
        N = pieceNum;
        headPVA = headState;
        tailPVA = tailState;
        T1.resize(N);
        A.create(6 * N, 6, 6);
        b.resize(6 * N, 3);
        gdC.resize(6 * N, 3);
        return;
    }

    inline void generate(const Eigen::MatrixXd &inPs,
                         const Eigen::VectorXd &ts)
    {
        T1 = ts;
        T2 = T1.cwiseProduct(T1);
        T3 = T2.cwiseProduct(T1);
        T4 = T2.cwiseProduct(T2);
        T5 = T4.cwiseProduct(T1);

        A.reset();
        b.setZero();

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        b.row(0) = headPVA.col(0).transpose();
        b.row(1) = headPVA.col(1).transpose();
        b.row(2) = headPVA.col(2).transpose();

        for (int i = 0; i < N - 1; i++)
        {
            A(6 * i + 3, 6 * i + 3) = 6.0;
            A(6 * i + 3, 6 * i + 4) = 24.0 * T1(i);
            A(6 * i + 3, 6 * i + 5) = 60.0 * T2(i);
            A(6 * i + 3, 6 * i + 9) = -6.0;
            A(6 * i + 4, 6 * i + 4) = 24.0;
            A(6 * i + 4, 6 * i + 5) = 120.0 * T1(i);
            A(6 * i + 4, 6 * i + 10) = -24.0;
            A(6 * i + 5, 6 * i) = 1.0;
            A(6 * i + 5, 6 * i + 1) = T1(i);
            A(6 * i + 5, 6 * i + 2) = T2(i);
            A(6 * i + 5, 6 * i + 3) = T3(i);
            A(6 * i + 5, 6 * i + 4) = T4(i);
            A(6 * i + 5, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i) = 1.0;
            A(6 * i + 6, 6 * i + 1) = T1(i);
            A(6 * i + 6, 6 * i + 2) = T2(i);
            A(6 * i + 6, 6 * i + 3) = T3(i);
            A(6 * i + 6, 6 * i + 4) = T4(i);
            A(6 * i + 6, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i + 6) = -1.0;
            A(6 * i + 7, 6 * i + 1) = 1.0;
            A(6 * i + 7, 6 * i + 2) = 2 * T1(i);
            A(6 * i + 7, 6 * i + 3) = 3 * T2(i);
            A(6 * i + 7, 6 * i + 4) = 4 * T3(i);
            A(6 * i + 7, 6 * i + 5) = 5 * T4(i);
            A(6 * i + 7, 6 * i + 7) = -1.0;
            A(6 * i + 8, 6 * i + 2) = 2.0;
            A(6 * i + 8, 6 * i + 3) = 6 * T1(i);
            A(6 * i + 8, 6 * i + 4) = 12 * T2(i);
            A(6 * i + 8, 6 * i + 5) = 20 * T3(i);
            A(6 * i + 8, 6 * i + 8) = -2.0;

            b.row(6 * i + 5) = inPs.col(i).transpose();
        }

        A(6 * N - 3, 6 * N - 6) = 1.0;
        A(6 * N - 3, 6 * N - 5) = T1(N - 1);
        A(6 * N - 3, 6 * N - 4) = T2(N - 1);
        A(6 * N - 3, 6 * N - 3) = T3(N - 1);
        A(6 * N - 3, 6 * N - 2) = T4(N - 1);
        A(6 * N - 3, 6 * N - 1) = T5(N - 1);
        A(6 * N - 2, 6 * N - 5) = 1.0;
        A(6 * N - 2, 6 * N - 4) = 2 * T1(N - 1);
        A(6 * N - 2, 6 * N - 3) = 3 * T2(N - 1);
        A(6 * N - 2, 6 * N - 2) = 4 * T3(N - 1);
        A(6 * N - 2, 6 * N - 1) = 5 * T4(N - 1);
        A(6 * N - 1, 6 * N - 4) = 2;
        A(6 * N - 1, 6 * N - 3) = 6 * T1(N - 1);
        A(6 * N - 1, 6 * N - 2) = 12 * T2(N - 1);
        A(6 * N - 1, 6 * N - 1) = 20 * T3(N - 1);

        b.row(6 * N - 3) = tailPVA.col(0).transpose();
        b.row(6 * N - 2) = tailPVA.col(1).transpose();
        b.row(6 * N - 1) = tailPVA.col(2).transpose();

        A.factorizeLU();
        A.solve(b);

        return;
    }

    inline double getTrajJerkCost() const
    {
        double objective = 0.0;
        for (int i = 0; i < N; i++)
        {
            objective += 36.0 * b.row(6 * i + 3).squaredNorm() * T1(i) +
                         144.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T2(i) +
                         192.0 * b.row(6 * i + 4).squaredNorm() * T3(i) +
                         240.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T3(i) +
                         720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T4(i) +
                         720.0 * b.row(6 * i + 5).squaredNorm() * T5(i);
        }
        return objective;
    }

    template <typename EIGENVEC, typename EIGENMAT>
    inline void evalTrajCostGrad(const Eigen::VectorXi &cons,
                                 const Eigen::VectorXi &idxHs,
                                 const std::vector<Eigen::MatrixXd> &cfgHs,
                                 const double &vmax,
                                 const double &amax,
                                 const Eigen::Vector3d &ci,
                                 double &cost,
                                 EIGENVEC &gdT,
                                 EIGENMAT &gdInPs)
    {
        gdT.setZero();
        gdInPs.setZero();
        gdC.setZero();

        cost = getTrajJerkCost();
        addGradJbyT(gdT);
        addGradJbyC(gdC);

        addTimeIntPenalty(cons, idxHs, cfgHs, vmax, amax, ci, cost, gdT, gdC);

        solveAdjGradC(gdC);
        addPropCtoT(gdC, gdT);
        addPropCtoP(gdC, gdInPs);
    }

    inline Trajectory getTraj(void) const
    {
        Trajectory traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++)
        {
            traj.emplace_back(T1(i), b.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
        }
        return traj;
    }


    // GaaiLam
    template <typename EIGENVEC, typename EIGENMAT>
    inline void initGradCost(EIGENVEC &gdT,
                             EIGENMAT &gdInPs, 
                             double &cost) {
        gdT.setZero();
        gdInPs.setZero();
        gdC.setZero();
        cost = getTrajJerkCost();
        addGradJbyT(gdT);
        addGradJbyC(gdC);
    }

    template <class TRAJGEN, typename EIGENVEC>
    // TRAJGEN::grad_cost_p(const Eigen::Vector3d &p, 
    //                      Eigen::Vector3d &gradp, 
    //                      double &cost) {}
    inline void addGrad2PVA(TRAJGEN *ptrObj,  
                            EIGENVEC &gdT, 
                            double &cost, 
                            const int &K=4) {
        //
        Eigen::Vector3d pos, vel, acc, jer;
        Eigen::Vector3d gradp, gradv, grada;
        double costp, costv, costa;
        Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3;
        double s1, s2, s3, s4, s5;
        double step, alpha;
        Eigen::Matrix<double, 6, 3> gradViolaPc, gradViolaVc, gradViolaAc;
        double gradViolaPt, gradViolaVt, gradViolaAt;
        double omg;

        int innerLoop;
        for (int i=0; i<N; ++i) {
            const auto &c = b.block<6, 3>(i * 6, 0);
            step = T1(i) / K;
            s1 = 0.0;
            innerLoop = K+1;

            for (int j=0; j<innerLoop; ++j) {
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                beta0 << 1.0, s1, s2, s3, s4, s5;
                beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
                beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
                beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
                alpha = 1.0 / K * j;
                pos = c.transpose() * beta0;
                vel = c.transpose() * beta1;
                acc = c.transpose() * beta2;
                jer = c.transpose() * beta3;

                omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;

                if ( ptrObj->grad_cost_p(pos, gradp, costp) ) {
                    gradViolaPc = beta0 * gradp.transpose();
                    gradViolaPt = alpha * gradp.transpose() * vel;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * gradViolaPc;
                    gdT(i) += omg * (costp/K + step * gradViolaPt);
                    cost += omg * step * costp;
                }
                if ( ptrObj->grad_cost_v(vel, gradv, costv) ) {
                    gradViolaVc = beta1 * gradv.transpose();
                    gradViolaVt = alpha * gradv.transpose() * acc;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * gradViolaVc;
                    gdT(i) += omg * (costv/K + step * gradViolaVt);
                    cost += omg * step * costv;
                }
                if ( ptrObj->grad_cost_a(acc, grada, costa) ) {
                    gradViolaAc = beta2 * grada.transpose();
                    gradViolaAt = alpha * grada.transpose() * jer;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * gradViolaAc;
                    gdT(i) += omg * (costa/K + step * gradViolaAt);
                    cost += omg * step * costa;
                }

                s1 += step;
            }
        }
    }

    template <class TRAJGEN, typename EIGENVEC>
    // i: the ith piece of traj
    // w: percent of time of this piece
    // TRAJGEN::grad_cost_p_at
    inline void addGrad2P_at (TRAJGEN *ptrObj, 
                              const int &i, 
                              const double &w, 
                              EIGENVEC &gdT, 
                              double &cost) {
        const auto &c = b.block<6, 3>(i * 6, 0);
        Eigen::Vector3d pos, vel, gradp;
        double costp = 0;
        double s1 = w * T1(i);
        double s2 = s1 * s1;
        double s3 = s2 * s1;
        double s4 = s3 * s1;
        double s5 = s4 * s1;
        Eigen::Matrix<double, 6, 1> beta0, beta1;
        beta0 << 1.0, s1, s2, s3, s4, s5;
        beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
        pos = c.transpose() * beta0;
        vel = c.transpose() * beta1;
        if ( ptrObj->grad_cost_p_at(pos, gradp, costp) ) {
            gdC.block<6, 3>(i * 6, 0) += beta0 * gradp.transpose();
            gdT(i) += w * gradp.transpose() * vel;
            cost += costp;
        }
    }

    template <typename EIGENVEC, typename EIGENMAT>
    inline void getGrad2TP(EIGENVEC &gdT,
                           EIGENMAT &gdInPs) {
        solveAdjGradC(gdC);
        addPropCtoT(gdC, gdT);
        addPropCtoP(gdC, gdInPs);
    }
};

#endif
