#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <limits>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// Constants
const double tolerance = 5e-3;
const int maxiter = 500;
const double pconverg = 1e-3;
const double ctbar = 138.0;

// Global Data equivalent to Fortran's MODULE hdata
struct HData {
  int ncomp;
  int nfiles;
  int tnum;
  vector<string> filename;
  vector<vector<int>> m_partf; // 2 x ncomp
  vector<vector<int>> m_part;  // 2 x ncomp
  double betap;
  vector<double> nliq;
  vector<double> ngas;
  vector<double> avgnum;
  double Z, Zgas, Zliq;
  vector<double> nmid;
  vector<double> slope;
  vector<double> beta;
  vector<vector<double>> mu; // ncomp x nfiles
  vector<double> mup;
  vector<int> nentry;
  vector<double> weight;
  vector<double> oldweight;

  int minp_global;
  int maxp_global;
  vector<vector<double>> dens; // ncomp x (maxp - minp + 1)

  vector<double> mw; // molecular weight of each component

  double energy_gas, energy_liq;

  vector<vector<vector<int>>> n; // ncomp x nfiles x nentry
  vector<vector<double>> e;      // nfiles x nentry
} hd;

// Parses one data line of a histogram file: ncomp integer particle numbers
// followed by the energy.  The whole line must be consumed, so a truncated,
// non-numeric or non-integer field is reported instead of silently ending
// the read.  Returns false on a malformed line.
bool parse_entry(const string &line, int ncomp, vector<int> &nvals,
                 double &energy) {
  istringstream ss(line);
  string tok;
  size_t pos;
  try {
    for (int c = 0; c < ncomp; ++c) {
      if (!(ss >> tok))
        return false;
      long v = stol(tok, &pos);
      if (pos != tok.size())
        return false;
      nvals[c] = (int)v;
    }
    if (!(ss >> tok))
      return false;
    energy = stod(tok, &pos);
    if (pos != tok.size())
      return false;
  } catch (const exception &) {
    return false;
  }
  return !(ss >> tok); // nothing may follow the energy
}

void Zcalc(double &nbelow) {
  nbelow = 0.0;
  hd.Z = 0.0;
  hd.Zgas = 0.0;
  hd.Zliq = 0.0;

  for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
    fill(hd.dens[icomp].begin(), hd.dens[icomp].end(), 0.0);
    hd.avgnum[icomp] = 0.0;
    hd.ngas[icomp] = 0.0;
    hd.nliq[icomp] = 0.0;
  }

  hd.energy_gas = 0.0;
  hd.energy_liq = 0.0;
  double split = 0.0;

  for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
    for (int i = 0; i < hd.nentry[ifile]; ++i) {
      double ylog = -1e9;
      for (int jfile = 0; jfile < hd.nfiles; ++jfile) {
        double prob = -(hd.beta[jfile] - hd.betap) * hd.e[ifile][i];
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          prob += (hd.beta[jfile] * hd.mu[icomp][jfile] -
                   hd.betap * hd.mup[icomp]) *
                  hd.n[icomp][ifile][i];
        }
        double tmp1 = prob + log(hd.nentry[jfile] / hd.oldweight[jfile]);
        ylog = max(ylog, tmp1) + log(1.0 + exp(-abs(ylog - tmp1)));
      }

      double exp_neg_ylog = exp(-ylog);
      hd.Z += exp_neg_ylog;
      split = 0.0;
      for (int icomp = 1; icomp < hd.ncomp; ++icomp) {
        split += hd.n[icomp][ifile][i] * hd.slope[hd.tnum];
      }

      if (hd.n[0][ifile][i] < hd.nmid[hd.tnum] + split) {
        hd.Zgas += exp_neg_ylog;
        nbelow += exp_neg_ylog;
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          hd.ngas[icomp] += hd.n[icomp][ifile][i] * exp_neg_ylog;
        }
        hd.energy_gas += hd.e[ifile][i] * exp_neg_ylog;
      } else {
        hd.Zliq += exp_neg_ylog;
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          hd.nliq[icomp] += hd.n[icomp][ifile][i] * exp_neg_ylog;
        }
        hd.energy_liq += hd.e[ifile][i] * exp_neg_ylog;
      }

      for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
        hd.avgnum[icomp] += hd.n[icomp][ifile][i] * exp_neg_ylog;
        int n_val = hd.n[icomp][ifile][i];
        int dens_idx = n_val - hd.minp_global;
        if (dens_idx >= 0 && dens_idx < hd.dens[icomp].size()) {
          hd.dens[icomp][dens_idx] += exp_neg_ylog;
        }
      }
    }
  }

  for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
    hd.ngas[icomp] /= hd.Z;
    hd.nliq[icomp] /= hd.Z;
    hd.avgnum[icomp] /= hd.Z;
    for (size_t i = 0; i < hd.dens[icomp].size(); ++i) {
      hd.dens[icomp][i] /= hd.Z;
    }
  }
  hd.energy_gas /= hd.Z;
  hd.energy_liq /= hd.Z;
  nbelow /= hd.Z;
}

// Gaussian elimination for Ax = b
// Solves in-place, modifying A and b. Returns false if singular.
bool gauss_solve(vector<vector<double>> &A, vector<double> &b,
                 vector<double> &x) {
  int n = A.size();
  for (int i = 0; i < n; i++) {
    double maxEl = fabs(A[i][i]);
    int maxRow = i;
    for (int k = i + 1; k < n; k++) {
      if (fabs(A[k][i]) > maxEl) {
        maxEl = fabs(A[k][i]);
        maxRow = k;
      }
    }
    if (maxEl < 1e-20)
      return false;

    for (int k = i; k < n; k++)
      swap(A[maxRow][k], A[i][k]);
    swap(b[maxRow], b[i]);

    for (int k = i + 1; k < n; k++) {
      double c = -A[k][i] / A[i][i];
      for (int j = i; j < n; j++) {
        if (i == j)
          A[k][j] = 0;
        else
          A[k][j] += c * A[i][j];
      }
      b[k] += c * b[i];
    }
  }
  x.assign(n, 0.0);
  for (int i = n - 1; i >= 0; i--) {
    x[i] = b[i];
    for (int k = i + 1; k < n; k++)
      x[i] -= A[i][k] * x[k];
    x[i] = x[i] / A[i][i];
  }
  return true;
}

// Refine hd.oldweight (read from input_fsp2.dat) by Newton-Raphson in log
// space; same solver as patch.cpp.  With g_j = ln(weight_j), the weights are
// the minimizer of the convex function (Shirts & Chodera, J. Chem. Phys. 129,
// 124105 (2008))
//   Phi(g) = sum_d log sum_j N_j exp(A_dj - g_j) + sum_j N_j g_j,
//   A_dj = -beta_j E_d + sum_c beta_j mu_cj N_cd,
// whose gradient is grad_k = N_k - W_k, with W_k / N_k = (fixed-point
// update) / (current weight).  Converged when maxd = max_k |W_k / N_k - 1|
// <= tol and the Newton step changes every ln(weight) by less than step_tol.
// On return hd.weight = hd.oldweight = converged weights, weight[0] = 1.
// Returns the number of iterations.
int solve_weights(double tol) {
  const int nf = hd.nfiles;
  const int max_iter = 1000;
  const double max_step = 1000.0; // largest allowed change in any ln(weight)
  const double step_tol = 1e-6;   // converged ln(weight) changes less than this

  long N_total = 0;
  for (int j = 0; j < nf; ++j)
    N_total += hd.nentry[j];

  vector<double> A(N_total * nf);
  long data_idx = 0;
  for (int ifile = 0; ifile < nf; ++ifile) {
    for (int i = 0; i < hd.nentry[ifile]; ++i) {
      for (int jfile = 0; jfile < nf; ++jfile) {
        double a_val = -hd.beta[jfile] * hd.e[ifile][i];
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          a_val += hd.beta[jfile] * hd.mu[icomp][jfile] * hd.n[icomp][ifile][i];
        }
        A[data_idx * nf + jfile] = a_val;
      }
      data_idx++;
    }
  }

  vector<double> g(nf), log_N(nf);
  for (int j = 0; j < nf; ++j) {
    if (!(hd.oldweight[j] > 0.0) || !isfinite(hd.oldweight[j])) {
      cerr << "Invalid weight " << hd.oldweight[j] << " for histogram "
           << hd.filename[j] << " in input_fsp2.dat\n";
      exit(1);
    }
    g[j] = log(hd.oldweight[j]) - log(hd.oldweight[0]);
    log_N[j] = log((double)hd.nentry[j]);
  }
  vector<double> denom_log(N_total);

  // Evaluates Phi(gg) and fills denom_log[d] = log sum_j N_j exp(A_dj - gg_j)
  auto objective = [&](const vector<double> &gg) {
    double phi = 0.0;
#pragma omp parallel for reduction(+ : phi)
    for (long d = 0; d < N_total; ++d) {
      double ymax = -1e300;
      for (int j = 0; j < nf; ++j)
        ymax = max(ymax, A[d * nf + j] + log_N[j] - gg[j]);
      double sum = 0.0;
      for (int j = 0; j < nf; ++j)
        sum += exp(A[d * nf + j] + log_N[j] - gg[j] - ymax);
      denom_log[d] = ymax + log(sum);
      phi += denom_log[d];
    }
    for (int j = 0; j < nf; ++j)
      phi += hd.nentry[j] * gg[j];
    return phi;
  };

  int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }
#endif

  double phi = objective(g);
  vector<double> W(nf), grad(nf), H(nf * nf);
  int iter = 0;
  double maxd = 0.0;

  while (true) {
    ++iter;

    // W_k = sum_d p_dk and Hessian H_kl = sum_d (p_dk delta_kl - p_dk p_dl),
    // with p_dk = N_k exp(A_dk - g_k) / denom_d.  Per-thread accumulation
    // followed by a deterministic reduction.
    vector<vector<double>> thread_W(nthreads, vector<double>(nf, 0.0));
    vector<vector<double>> thread_H(nthreads, vector<double>(nf * nf, 0.0));
#pragma omp parallel
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      vector<double> lW(nf, 0.0), lH(nf * nf, 0.0), p(nf);
#pragma omp for schedule(static)
      for (long d = 0; d < N_total; ++d) {
        for (int k = 0; k < nf; ++k) {
          p[k] = exp(A[d * nf + k] + log_N[k] - g[k] - denom_log[d]);
          lW[k] += p[k];
        }
        // Diagonal written as p_k * sum_{j!=k} p_j rather than p_k - p_k^2,
        // which cancels catastrophically when one histogram has p_k ~ 1
        for (int k = 0; k < nf; ++k)
          for (int l = 0; l < nf; ++l)
            if (l != k) {
              lH[k * nf + l] -= p[k] * p[l];
              lH[k * nf + k] += p[k] * p[l];
            }
      }
      thread_W[tid] = lW;
      thread_H[tid] = lH;
    }
    fill(W.begin(), W.end(), 0.0);
    fill(H.begin(), H.end(), 0.0);
    for (int t = 0; t < nthreads; ++t) {
      for (int k = 0; k < nf; ++k)
        W[k] += thread_W[t][k];
      for (int kl = 0; kl < nf * nf; ++kl)
        H[kl] += thread_H[t][kl];
    }
    maxd = 0.0;
    for (int k = 0; k < nf; ++k) {
      grad[k] = hd.nentry[k] - W[k];
      maxd = max(maxd, fabs(W[k] / hd.nentry[k] - 1.0));
    }
    if (!isfinite(maxd)) {
      cerr << "Error: non-finite weights at iteration " << iter << endl;
      exit(1);
    }

    // Newton step on the reduced system (g_0 held fixed)
    int nr = nf - 1;
    vector<vector<double>> Hr(nr, vector<double>(nr));
    vector<double> rhs(nr), step_r;
    for (int k = 0; k < nr; ++k) {
      rhs[k] = -grad[k + 1];
      for (int l = 0; l < nr; ++l)
        Hr[k][l] = H[(k + 1) * nf + (l + 1)];
    }
    vector<double> step(nf, 0.0);
    bool newton_ok = gauss_solve(Hr, rhs, step_r);
    double max_newton = 0.0;
    if (newton_ok) {
      for (int k = 0; k < nr; ++k) {
        step[k + 1] = step_r[k];
        max_newton = max(max_newton, fabs(step_r[k]));
      }
    }

    cout << "iteration = " << setw(5) << iter << "     deviation = "
         << scientific << setprecision(4) << maxd;
    if (newton_ok)
      cout << "     Newton step = " << max_newton;
    cout << fixed << endl;

    // A small maxd alone is not enough: a weakly overlapping histogram can
    // still be off in ln(weight), so also require the Newton step to vanish
    if ((maxd <= tol && newton_ok && max_newton <= step_tol) ||
        iter >= max_iter)
      break;

    // Backtracking line search on Phi along step (Armijo condition).
    // Updates g, phi and denom_log and returns true on success.
    vector<double> g_new(nf);
    auto line_search = [&](vector<double> step) {
      double big = 0.0;
      for (int k = 0; k < nf; ++k)
        big = max(big, fabs(step[k]));
      if (big > max_step)
        for (int k = 0; k < nf; ++k)
          step[k] *= max_step / big;
      double slope = 0.0;
      for (int k = 0; k < nf; ++k)
        slope += grad[k] * step[k];
      if (!(slope < 0.0))
        return false; // not a descent direction
      double lambda = 1.0;
      for (int ls = 0; ls < 40; ++ls) {
        for (int k = 0; k < nf; ++k)
          g_new[k] = g[k] + lambda * step[k];
        double phi_new = objective(g_new);
        if (phi_new <= phi + 1e-4 * lambda * slope + 1e-12 * fabs(phi)) {
          g = g_new;
          phi = phi_new; // denom_log now corresponds to g
          return true;
        }
        lambda *= 0.5;
      }
      objective(g); // restore denom_log for the unchanged g
      return false;
    };

    bool moved = false;
    if (newton_ok)
      moved = line_search(step);
    if (!moved) {
      // Newton failed (nearly singular Hessian far from the solution, where
      // histograms do not yet overlap): take a fixed-point step instead
      for (int k = 1; k < nf; ++k)
        step[k] = log(max(W[k], 1e-300) / hd.nentry[k]) -
                  log(max(W[0], 1e-300) / hd.nentry[0]);
      moved = line_search(step);
    }
    if (!moved) {
      cerr << "Error: line search failed at iteration " << iter << endl;
      break;
    }
  }

  if (iter >= max_iter) {
    cerr << "Warning: weights did not converge in " << max_iter
         << " iterations (deviation = " << scientific << maxd << ")" << fixed
         << endl;
  }
  for (int j = 0; j < nf; ++j) {
    hd.weight[j] = exp(g[j]);
    hd.oldweight[j] = hd.weight[j];
  }
  return iter;
}

int main() {
  bool lphase = true;
  bool lpvt = false;
  bool lfixp = false;

  cout << "PHASE VERSION 8.0" << endl;
  if (lphase && lpvt) {
    cout << "Can not do phase coexistence and PVT at the same time\nProgram "
            "terminating\n";
    return 1;
  } else if (!lphase && !lpvt) {
    cout << "Must set lphase or lpvt = true\nProgram terminating\n";
    return 1;
  }

  if (lphase) {
    cout << "PHASE COEXISTENCE MODE\n";
    if (lfixp)
      cout << "\nFIXED PRESSURE MODE\n";
  } else if (lpvt) {
    cout << "PVT MODE\n";
  }

  hd.nmid.clear();
  hd.slope.clear();

  ifstream f_input("input_fsp2.dat");
  if (!f_input) {
    cerr << "Could not open input_fsp2.dat\n";
    return 1;
  }
  string suffix;
  f_input >> hd.ncomp >> hd.nfiles >> suffix;
  if (!f_input || hd.ncomp < 1 || hd.nfiles < 1) {
    cerr << "input_fsp2.dat must start with the number of components, the "
            "number of histograms and the file suffix\n";
    return 1;
  }

  vector<double> dmup(hd.ncomp);
  hd.avgnum.resize(hd.ncomp);
  hd.ngas.resize(hd.ncomp);
  hd.nliq.resize(hd.ncomp);

  hd.m_partf.assign(2, vector<int>(hd.ncomp, 0));
  hd.m_part.assign(2, vector<int>(hd.ncomp, 0));
  for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
    hd.m_part[0][icomp] = numeric_limits<int>::max();
    hd.m_part[1][icomp] = -1;
  }

  hd.filename.resize(hd.nfiles);
  hd.nentry.resize(hd.nfiles);
  hd.oldweight.resize(hd.nfiles);
  hd.beta.resize(hd.nfiles);
  hd.mu.assign(hd.ncomp, vector<double>(hd.nfiles, 0.0));
  hd.mup.resize(hd.ncomp);
  vector<vector<double>> mu_temp;
  hd.weight.resize(hd.nfiles);

  for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
    f_input >> hd.filename[ifile] >> hd.nentry[ifile] >> hd.oldweight[ifile] >>
        hd.beta[ifile];
    for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
      f_input >> hd.mu[icomp][ifile];
    }
    if (!f_input) {
      cerr << "input_fsp2.dat specifies " << hd.nfiles
           << " histograms but entry " << ifile + 1
           << " is missing or malformed\n";
      return 1;
    }
    hd.beta[ifile] = 1.0 / hd.beta[ifile];
  }
  f_input.close();

  hd.n.assign(hd.ncomp, vector<vector<int>>(hd.nfiles));
  hd.e.assign(hd.nfiles, vector<double>());

  double eng_min = 1e5;
  double eng_max = -1e5;

  cout << "histogram     npoints  tstar  chemical potentials  min max n1 n2 "
          "etc\n";
  cout << "===================================================================="
          "========\n";

  double xdim = 0, ydim = 0, zdim = 0;
  for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
    for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
      hd.m_partf[0][icomp] = numeric_limits<int>::max();
      hd.m_partf[1][icomp] = -1;
    }

    string hname = "his" + hd.filename[ifile] + suffix + ".dat";
    ifstream f_his(hname);
    if (!f_his) {
      cerr << "Could not open " << hname << endl;
      return 1;
    }

    // header: T, number of components, chemical potentials and box
    // dimensions, all on the first line
    double beta_val = 0.0;
    int ncomp2 = 0;
    double xdim1 = 0, ydim1 = 0, zdim1 = 0;
    string line;
    getline(f_his, line);
    istringstream header(line);
    header >> beta_val >> ncomp2;
    if (header && ncomp2 == hd.ncomp) {
      double mu_tmp;
      for (int icomp = 0; icomp < ncomp2; ++icomp)
        header >> mu_tmp;
      header >> xdim1 >> ydim1 >> zdim1;
    }
    if (!header) {
      cerr << "Could not read the header (line 1) of " << hname << endl;
      return 1;
    }

    if (hd.ncomp != ncomp2) {
      cout << "Number of components do not match in input_hs.dat and "
              "histogram\n";
      return 1;
    }

    hd.beta[ifile] = 1.0 / beta_val;
    if (ifile == 0) {
      xdim = xdim1;
      ydim = ydim1;
      zdim = zdim1;
    } else {
      if (xdim1 != xdim || ydim1 != ydim || zdim1 != zdim) {
        cout << "System size or temperature in file " << hname
             << " inconsistent with first file\n";
        return 1;
      }
    }

    vector<vector<int>> n_temp(hd.ncomp);
    vector<double> e_temp;

    vector<int> nt(hd.ncomp);
    double et;
    int lineno = 1;
    while (getline(f_his, line)) {
      ++lineno;
      if (line.find_first_not_of(" \t\r") == string::npos)
        continue; // blank line
      if (!parse_entry(line, hd.ncomp, nt, et)) {
        cerr << "Malformed entry on line " << lineno << " of " << hname
             << ":\n  \"" << line << "\"\nExpected " << hd.ncomp
             << " integer particle number(s) followed by the energy\n";
        return 1;
      }
      if (f_his.eof()) {
        // getline hit end of file without a newline: GOMC terminates every
        // line, so the file may have been cut off mid-write
        cerr << "Warning: last line of " << hname
             << " has no newline; the file may be truncated\n";
      }

      e_temp.push_back(et);
      for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
        n_temp[icomp].push_back(nt[icomp]);
        hd.m_part[0][icomp] = min(hd.m_part[0][icomp], nt[icomp]);
        hd.m_part[1][icomp] = max(hd.m_part[1][icomp], nt[icomp]);
        hd.m_partf[0][icomp] = min(hd.m_partf[0][icomp], nt[icomp]);
        hd.m_partf[1][icomp] = max(hd.m_partf[1][icomp], nt[icomp]);
      }
      eng_min = min(eng_min, et);
      eng_max = max(eng_max, et);
    }

    if (e_temp.empty()) {
      cerr << "Histogram file " << hname << " contains no entries\n";
      return 1;
    }
    hd.nentry[ifile] = e_temp.size();
    for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
      hd.n[icomp][ifile] = n_temp[icomp];
    }
    hd.e[ifile] = e_temp;

    cout << hname << " " << hd.nentry[ifile] << " " << fixed << setprecision(2)
         << 1.0 / hd.beta[ifile] << " ";
    for (int icomp = 0; icomp < hd.ncomp; ++icomp)
      cout << hd.mu[icomp][ifile] << " ";
    for (int icomp = 0; icomp < hd.ncomp; ++icomp)
      cout << hd.m_partf[0][icomp] << " " << hd.m_partf[1][icomp] << " ";
    cout << endl;
  }

  dmup[0] = hd.mu[0][0] / 10000.0;
  double vol = xdim * ydim * zdim;

  if (lfixp) {
    if (hd.ncomp < 2) {
      cout << "FIXED PRESSURE REQUIRES NCOMP > 1\nPROGRAM TERMINATING\n";
      return 1;
    }
    dmup[hd.ncomp - 1] = hd.mu[hd.ncomp - 1][0] / 1000.0;
  }

  hd.minp_global = 100000;
  hd.maxp_global = -100000;
  for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
    hd.minp_global = min(hd.m_part[0][icomp], hd.minp_global);
    hd.maxp_global = max(hd.m_part[1][icomp], hd.maxp_global);
  }

  hd.dens.assign(hd.ncomp,
                 vector<double>(hd.maxp_global - hd.minp_global + 1, 0.0));

  int iter = solve_weights(tolerance);

  ofstream f_out24("input_fsp2.dat");
  f_out24 << hd.ncomp << "\n" << hd.nfiles << "\n" << suffix << "\n";
  for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
    f_out24 << setw(5) << hd.filename[ifile] << setw(10) << hd.nentry[ifile]
            << " " << scientific << setprecision(6) << hd.weight[ifile] << " "
            << fixed << setprecision(4) << 1.0 / hd.beta[ifile] << " ";
    for (int icomp = 0; icomp < hd.ncomp; ++icomp)
      f_out24 << hd.mu[icomp][ifile] << " ";
    f_out24 << "\n";
  }
  f_out24 << "total iterations = " << iter << "\n";
  f_out24.close();

  vector<double> t_new, mu2min, mu2max, mu2_incr;
  double lnZ0 = 0.0, pset = 0.0;
  int ntemp = 0;

  if (lphase) {
    cout << "Starting coexistence calculation\n";
    ofstream f_phase("phase.dat");
    f_phase << "# Suffixes and Files = " << suffix;
    for (int ifile = 0; ifile < hd.nfiles; ++ifile)
      f_phase << hd.filename[ifile] << " ";
    f_phase << "\n#" << setw(12) << "T";
    for (int icomp = 0; icomp < hd.ncomp; ++icomp)
      f_phase << setw(14) << "mu" + to_string(icomp + 1);
    f_phase << setw(14) << "Eng_liq" << setw(14) << "Eng_gas";
    for (int icomp = 0; icomp < hd.ncomp; ++icomp)
      f_phase << setw(14) << "N_liq" + to_string(icomp + 1);
    for (int icomp = 0; icomp < hd.ncomp; ++icomp)
      f_phase << setw(14) << "N_gas" + to_string(icomp + 1);
    f_phase << setw(14) << "Rho_liq" << setw(14) << "Rho_gas";
    f_phase << setw(14) << "ln(Zliq)" << setw(14) << "ln(Zgas)\n";
    f_phase.close();

    ofstream f_trho_liq("trho_liq.dat");
    f_trho_liq << "# Density(kg/m^3)   Tstar\n";
    f_trho_liq.close();

    ofstream f_trho_vap("trho_vap.dat");
    f_trho_vap << "# Density(kg/m^3)   Tstar\n";
    f_trho_vap.close();

    ifstream f_ph("phinput.idat");
    if (!f_ph.is_open()) {
      cout << "Warning: Could not open phinput.idat\n";
    } else {
      cout << "Successfully opened phinput.idat. Here are its contents:\n";
      cout << "--------------------------------------------------\n";
      string line;
      while (getline(f_ph, line)) {
        cout << line << endl;
      }
      cout << "--------------------------------------------------\n";

      // Re-open or clear the stream to parse properly
      f_ph.clear();
      f_ph.seekg(0);

      // Line 1: Header
      getline(f_ph, line);

      // Line 2: Initial guess variables
      getline(f_ph, line);
      stringstream ss_mu(line);
      if (lfixp) {
        if (hd.ncomp > 2) {
          ss_mu >> hd.mup[0] >> hd.mup[hd.ncomp - 2] >> lnZ0 >> pset;
        } else {
          ss_mu >> hd.mup[0] >> hd.mup[hd.ncomp - 1] >> lnZ0 >> pset;
        }
      } else {
        ss_mu >> hd.mup[0] >> lnZ0;
      }

      // Lines 3+: Read MW if present, and any valid 6-parameter conditions
      hd.mw.assign(hd.ncomp, 1.0); // default MW if not provided
      bool next_line_is_mw = false;

      while (getline(f_ph, line)) {
        if (next_line_is_mw) {
          stringstream ss_mw(line);
          for (int i = 0; i < hd.ncomp; ++i) {
            ss_mw >> hd.mw[i];
          }
          next_line_is_mw = false;
          continue;
        }

        if (line.find("MW") != string::npos) {
          next_line_is_mw = true;
          continue;
        }

        stringstream ss(line);
        double t, m2min, m2max, m2inc, nmid_val, slope_val;
        if (ss >> t >> m2min >> m2max >> m2inc >> nmid_val >> slope_val) {
          t_new.push_back(t);
          mu2min.push_back(m2min);
          mu2max.push_back(m2max);
          mu2_incr.push_back(m2inc);
          hd.nmid.push_back(nmid_val);
          hd.slope.push_back(slope_val);
        }
      }
      ntemp = t_new.size();
      f_ph.close();

      cout << "Read " << ntemp << " conditions from phinput.idat:" << endl;
      cout << "  lnZ0=" << lnZ0 << endl;
      if (hd.ncomp > 0) {
        cout << "  Parsed Molecular Weights: ";
        for (int i = 0; i < hd.ncomp; ++i)
          cout << hd.mw[i] << " ";
        cout << endl;
      }
      for (int i = 0; i < ntemp; ++i) {
        cout << "  T=" << t_new[i] << " m2min=" << mu2min[i]
             << " m2max=" << mu2max[i] << " m2inc=" << mu2_incr[i]
             << " nmid=" << hd.nmid[i] << " slope=" << hd.slope[i] << endl;
      }
    }

  } else if (lpvt) {
    ofstream f_pvt("pvt.dat");
    f_pvt << "# Suffixes and Files = " << suffix;
    for (int ifile = 0; ifile < hd.nfiles; ++ifile)
      f_pvt << hd.filename[ifile] << " ";
    f_pvt << "\n#    mu1          mu2           <N>           lnZ\n";
    f_pvt.close();

    ifstream f_pvtin("pvt.idat");
    double t;
    while (f_pvtin >> t) {
      t_new.push_back(t);
      vector<double> m_tmp(hd.ncomp);
      for (int icomp = 0; icomp < hd.ncomp; ++icomp)
        f_pvtin >> m_tmp[icomp];
      mu_temp.push_back(m_tmp);
    }
    ntemp = t_new.size();
    f_pvtin.close();
  }

  for (hd.tnum = 0; hd.tnum < ntemp; ++hd.tnum) {
    hd.betap = 1.0 / t_new[hd.tnum];
    int num = 0;
    iter = 0;

    if (lphase) {
      if (mu2_incr[hd.tnum] > 0 && hd.ncomp > 2)
        hd.mup[hd.ncomp - 1] = mu2min[hd.tnum];
      else if (mu2_incr[hd.tnum] < 0 && hd.ncomp > 2)
        hd.mup[hd.ncomp - 1] = mu2max[hd.tnum];
      else if (hd.ncomp == 2)
        hd.mup[hd.ncomp - 1] = mu2min[hd.tnum];

      bool keep_going = true;
      while (keep_going) {
        if (hd.ncomp == 1 && iter == 0)
          keep_going = true;
        else if (hd.mup[hd.ncomp - 1] < mu2max[hd.tnum] + 0.00001 &&
                 hd.mup[hd.ncomp - 1] > mu2min[hd.tnum] - 0.00001)
          keep_going = true;
        else
          keep_going = false;

        if (!keep_going)
          break;

        if (lfixp) {
          if (hd.ncomp < 2) {
            cout << "FIXED PRESSURE REQUIRES NCOMP > 1\n";
            return 1;
          }
          dmup[hd.ncomp - 1] = hd.mu[hd.ncomp - 1][0] / 1000.0;
        }

        int iter2 = 0;
        int flip2 = 0;
        double oldsign2 = 1.0;
        double diff2 = 1e6;
        num++;
        double press = 0.0;
        double nbelow = 0.0;
        double nbelow1 = 0.0;

        while (diff2 > pconverg && iter2 < maxiter) {
          iter = 0;

          while (abs(0.5 - nbelow) > 1e-6 && iter <= maxiter) {
            if (iter > 0) {
              hd.mup[0] += dmup[0];
              Zcalc(nbelow1);
              hd.mup[0] -= (log(nbelow1) - log(1.0 - nbelow1)) * dmup[0] /
                           (log(nbelow1) - log(nbelow) - log(1.0 - nbelow1) +
                            log(1.0 - nbelow));
            }
            Zcalc(nbelow);
            iter++;
          }

          press = (log(hd.Zliq) - lnZ0) * ctbar / hd.betap / vol;

          if (!lfixp) {
            diff2 = 1e-6;
          } else {
            double sign2 = press - pset;
            diff2 = abs(press - pset);

            if (flip2 > 0) {
              dmup[hd.ncomp - 1] /= 2.0;
              flip2 = 0;
            }

            double newsign2;
            if (sign2 > 0) {
              hd.mup[hd.ncomp - 1] += dmup[hd.ncomp - 1];
              newsign2 = -1;
            } else {
              hd.mup[hd.ncomp - 1] -= dmup[hd.ncomp - 1];
              newsign2 = 1;
            }

            if (newsign2 != oldsign2) {
              oldsign2 = newsign2;
              flip2++;
            }
          }
          iter2++;
        }

        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          hd.ngas[icomp] /= nbelow;
          hd.nliq[icomp] /= (1.0 - nbelow);
        }
        hd.energy_gas /= nbelow;
        hd.energy_liq /= (1.0 - nbelow);

        double rho_liq = 0.0, rho_gas = 0.0;
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          rho_liq += hd.nliq[icomp] * hd.mw[icomp];
          rho_gas += hd.ngas[icomp] * hd.mw[icomp];
        }
        // Convert N * (g/mol) / A^3 to kg/m^3
        rho_liq = (rho_liq / vol) * 1660.53906717;
        rho_gas = (rho_gas / vol) * 1660.53906717;

        ofstream f_phase("phase.dat", ios::app);
        f_phase << fixed << setprecision(4);
        f_phase << "  " << setw(12) << 1.0 / hd.betap;
        for (int icomp = 0; icomp < hd.ncomp; ++icomp)
          f_phase << setw(14) << hd.mup[icomp];
        f_phase << setw(14) << hd.energy_liq << setw(14) << hd.energy_gas;
        for (int icomp = 0; icomp < hd.ncomp; ++icomp)
          f_phase << setw(14) << hd.nliq[icomp];
        for (int icomp = 0; icomp < hd.ncomp; ++icomp)
          f_phase << setw(14) << hd.ngas[icomp];
        f_phase << setw(14) << rho_liq << setw(14) << rho_gas;
        f_phase << setw(14) << log(hd.Zliq) << setw(14) << log(hd.Zgas) << "\n";
        f_phase.close();

        double larger_rho = std::max(rho_liq, rho_gas);
        double smaller_rho = std::min(rho_liq, rho_gas);

        ofstream f_tliq("trho_liq.dat", ios::app);
        f_tliq << fixed << setprecision(4) << setw(16) << larger_rho << setw(12)
               << 1.0 / hd.betap << "\n";
        f_tliq.close();

        ofstream f_tvap("trho_vap.dat", ios::app);
        f_tvap << fixed << setprecision(4) << setw(16) << smaller_rho
               << setw(12) << 1.0 / hd.betap << "\n";
        f_tvap.close();

        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
          string fname3 = "d" + to_string(icomp + 1) + "n" +
                          to_string(hd.tnum + 1) + "t" + to_string(num) +
                          "a.dat";
          ofstream f_d(fname3);
          f_d << "# " << fixed << setprecision(3) << 1.0 / hd.betap << " ";
          for (int jcomp = 0; jcomp < hd.ncomp; ++jcomp)
            f_d << fixed << setprecision(2) << hd.mup[jcomp] << " ";
          f_d << xdim << " " << ydim << " " << zdim << "\n";

          for (int ipart = hd.m_part[0][icomp]; ipart <= hd.m_part[1][icomp];
               ++ipart) {
            int dens_idx = ipart - hd.minp_global;
            if (dens_idx >= 0 && dens_idx < hd.dens[icomp].size()) {
              f_d << setw(6) << ipart << " " << scientific << setprecision(4)
                  << hd.dens[icomp][dens_idx] << "\n";
            }
          }
          f_d.close();
        }

        if (hd.ncomp >= 2) {
          hd.mup[1] += mu2_incr[hd.tnum];
        }
      } // while (keep_going)

    } else if (lpvt) {
      num = 1;
      for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
        hd.mup[icomp] = mu_temp[hd.tnum][icomp];
      }
      hd.nmid[hd.tnum] = hd.m_part[1][0]; // equivalent to m_part(2,1)
      double nbelow;
      Zcalc(nbelow);

      double avg = 0.0;
      for (int icomp = 0; icomp < hd.ncomp; ++icomp)
        avg += hd.avgnum[icomp];

      ofstream f_pvt("pvt.dat", ios::app);
      f_pvt << fixed << setprecision(6) << avg << "     " << log(hd.Z) << "\n";
      f_pvt.close();

      for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
        string fname3 = "d" + to_string(icomp + 1) + "n" +
                        to_string(hd.tnum + 1) + "t" + to_string(num) + "a.dat";
        ofstream f_d(fname3);
        f_d << "# " << fixed << setprecision(3) << 1.0 / hd.betap << " ";
        for (int jcomp = 0; jcomp < hd.ncomp; ++jcomp)
          f_d << fixed << setprecision(2) << hd.mup[jcomp] << " ";
        f_d << xdim << " " << ydim << " " << zdim << "\n";

        for (int ipart = hd.m_part[0][icomp]; ipart <= hd.m_part[1][icomp];
             ++ipart) {
          int dens_idx = ipart - hd.minp_global;
          if (dens_idx >= 0 && dens_idx < hd.dens[icomp].size()) {
            f_d << setw(6) << ipart << " " << scientific << setprecision(4)
                << hd.dens[icomp][dens_idx] << "\n";
          }
        }
        f_d.close();
      }
    }
  }

  cout << "Phase coexistence calculation complete\nPROGRAM TERMINATING\n";
  return 0;
}
