#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdlib>

#include <limits>

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
    
    double energy_gas, energy_liq;

    vector<vector<vector<int>>> n; // ncomp x nfiles x nentry
    vector<vector<double>> e;      // nfiles x nentry
} hd;

void Zcalc(double& nbelow) {
    nbelow = 0.0;
    hd.Z = 0.0;
    hd.Zgas = 0.0;
    hd.Zliq = 0.0;
    
    for(int icomp=0; icomp<hd.ncomp; ++icomp) {
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
                    prob += (hd.beta[jfile] * hd.mu[icomp][jfile] - hd.betap * hd.mup[icomp]) * hd.n[icomp][ifile][i];
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
                if(dens_idx >= 0 && dens_idx < hd.dens[icomp].size()) {
                    hd.dens[icomp][dens_idx] += exp_neg_ylog;
                }
            }
        }
    }

    for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
        hd.ngas[icomp] /= hd.Z;
        hd.nliq[icomp] /= hd.Z;
        hd.avgnum[icomp] /= hd.Z;
        for (size_t i=0; i<hd.dens[icomp].size(); ++i) {
            hd.dens[icomp][i] /= hd.Z;
        }
    }
    hd.energy_gas /= hd.Z;
    hd.energy_liq /= hd.Z;
    nbelow /= hd.Z;
}

int main() {
    bool lphase = true;
    bool lpvt = false;
    bool lfixp = false;

    cout << "PHASE VERSION 8.0" << endl;
    if (lphase && lpvt) {
        cout << "Can not do phase coexistence and PVT at the same time\nProgram terminating\n";
        return 1;
    } else if (!lphase && !lpvt) {
        cout << "Must set lphase or lpvt = true\nProgram terminating\n";
        return 1;
    }

    if (lphase) {
        cout << "PHASE COEXISTENCE MODE\n";
        if (lfixp) cout << "\nFIXED PRESSURE MODE\n";
    } else if (lpvt) {
        cout << "PVT MODE\n";
    }

    hd.nmid.clear();
    hd.slope.clear();

    ifstream f_input("input_fsp2.dat");
    if(!f_input) { cerr << "Could not open input_fsp2.dat\n"; return 1; }
    string suffix;
    f_input >> hd.ncomp >> hd.nfiles >> suffix;

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
        f_input >> hd.filename[ifile] >> hd.nentry[ifile] >> hd.oldweight[ifile] >> hd.beta[ifile];
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
            f_input >> hd.mu[icomp][ifile];
        }
        hd.beta[ifile] = 1.0 / hd.beta[ifile];
    }
    f_input.close();

    hd.n.assign(hd.ncomp, vector<vector<int>>(hd.nfiles));
    hd.e.assign(hd.nfiles, vector<double>());

    double eng_min = 1e5;
    double eng_max = -1e5;

    cout << "histogram     npoints  tstar  chemical potentials  min max n1 n2 etc\n";
    cout << "============================================================================\n";

    double xdim = 0, ydim = 0, zdim = 0;
    for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
            hd.m_partf[0][icomp] = numeric_limits<int>::max();
            hd.m_partf[1][icomp] = -1;
        }

        string hname = "his" + hd.filename[ifile] + suffix + ".dat";
        ifstream f_his(hname);
        if(!f_his) { cerr << "Could not open " << hname << endl; return 1; }

        double beta_val;
        int ncomp2;
        f_his >> beta_val >> ncomp2;
        vector<double> mu_tmp(ncomp2);
        for(int icomp=0; icomp<ncomp2; ++icomp) f_his >> mu_tmp[icomp];
        double xdim1, ydim1, zdim1;
        f_his >> xdim1 >> ydim1 >> zdim1;

        if (hd.ncomp != ncomp2) {
            cout << "Number of components do not match in input_hs.dat and histogram\n";
            return 1;
        }

        hd.beta[ifile] = 1.0 / beta_val;
        if (ifile == 0) {
            xdim = xdim1; ydim = ydim1; zdim = zdim1;
        } else {
            if (xdim1 != xdim || ydim1 != ydim || zdim1 != zdim) {
                cout << "System size or temperature in file " << hname << " inconsistent with first file\n";
                return 1;
            }
        }

        vector<vector<int>> n_temp(hd.ncomp);
        vector<double> e_temp;

        while (true) {
            vector<int> nt(hd.ncomp);
            bool ok = true;
            for(int icomp=0; icomp<hd.ncomp; ++icomp) {
                if(!(f_his >> nt[icomp])) { ok = false; break; }
            }
            double et;
            if(!ok || !(f_his >> et)) break;

            e_temp.push_back(et);
            for(int icomp=0; icomp<hd.ncomp; ++icomp) {
                n_temp[icomp].push_back(nt[icomp]);
                hd.m_part[0][icomp] = min(hd.m_part[0][icomp], nt[icomp]);
                hd.m_part[1][icomp] = max(hd.m_part[1][icomp], nt[icomp]);
                hd.m_partf[0][icomp] = min(hd.m_partf[0][icomp], nt[icomp]);
                hd.m_partf[1][icomp] = max(hd.m_partf[1][icomp], nt[icomp]);
            }
            eng_min = min(eng_min, et);
            eng_max = max(eng_max, et);
        }

        hd.nentry[ifile] = e_temp.size();
        for(int icomp=0; icomp<hd.ncomp; ++icomp) {
            hd.n[icomp][ifile] = n_temp[icomp];
        }
        hd.e[ifile] = e_temp;

        cout << hname << " " << hd.nentry[ifile] << " " << fixed << setprecision(2) << 1.0/hd.beta[ifile] << " ";
        for(int icomp=0; icomp<hd.ncomp; ++icomp) cout << hd.mu[icomp][ifile] << " ";
        for(int icomp=0; icomp<hd.ncomp; ++icomp) cout << hd.m_partf[0][icomp] << " " << hd.m_partf[1][icomp] << " ";
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

    hd.dens.assign(hd.ncomp, vector<double>(hd.maxp_global - hd.minp_global + 1, 0.0));

    int iter = 0;
    double maxd = 100.0;
    while (maxd > tolerance && iter < 1e4) {
        iter++;
        maxd = 0.0;
        for (int kfile = 0; kfile < hd.nfiles; ++kfile) {
            hd.weight[kfile] = 0.0;
            for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
                for (int i = 0; i < hd.nentry[ifile]; ++i) {
                    double ylog = -1e9;
                    for (int jfile = 0; jfile < hd.nfiles; ++jfile) {
                        double prob = -(hd.beta[jfile] - hd.beta[kfile]) * hd.e[ifile][i];
                        for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
                            prob += (hd.beta[jfile] * hd.mu[icomp][jfile] - hd.beta[kfile] * hd.mu[icomp][kfile]) * hd.n[icomp][ifile][i];
                        }
                        double tmp1 = prob + log(hd.nentry[jfile] / hd.oldweight[jfile]);
                        ylog = max(ylog, tmp1) + log(1.0 + exp(-abs(ylog - tmp1)));
                    }
                    hd.weight[kfile] += exp(-ylog);
                }
            }
            double val = abs(hd.weight[kfile] / hd.oldweight[kfile] - 1.0);
            if (val > maxd) maxd = val;
        }
        for (int i = 1; i < hd.nfiles; ++i) {
            hd.weight[i] = hd.weight[i] / hd.weight[0];
        }
        hd.weight[0] = 1.0;
        for (int i = 0; i < hd.nfiles; ++i) {
            hd.oldweight[i] = hd.weight[i];
        }
        cout << "iteration = " << setw(5) << iter << "     deviation = " << fixed << setprecision(6) << maxd << endl;
    }

    ofstream f_out24("input_fsp2.dat");
    f_out24 << hd.ncomp << "\n" << hd.nfiles << "\n" << suffix << "\n";
    for (int ifile = 0; ifile < hd.nfiles; ++ifile) {
        f_out24 << setw(5) << hd.filename[ifile] << setw(10) << hd.nentry[ifile] << " " 
                << scientific << setprecision(6) << hd.weight[ifile] << " " 
                << fixed << setprecision(4) << 1.0/hd.beta[ifile] << " ";
        for(int icomp=0; icomp<hd.ncomp; ++icomp) f_out24 << hd.mu[icomp][ifile] << " ";
        f_out24 << "\n";
    }
    f_out24 << "total iterations = " << iter << "\n";
    f_out24.close();

    vector<double> t_new, mu2min, mu2max, mu2_incr;
    double lnZ0 = 0.0, pset = 0.0;
    int ntemp = 0;

    if (lphase) {
        cout << "Starting coexistence calculation\n";
        ofstream f_phase("phase.dat", ios::app);
        f_phase << "/* Suffixes and Files = " << suffix;
        for(int ifile=0; ifile<hd.nfiles; ++ifile) f_phase << hd.filename[ifile] << " ";
        f_phase << "\n/*      T            mu1           mu2        Eng_liq         Eng_gas       N_liq1          N_liq2        N_gas1         Ngas2  ln(Zliq) ln(Zgas)\n";
        f_phase.close();

        ifstream f_ph("phinput.idat");
        string dummy;
        if (lfixp) {
            getline(f_ph, dummy);
            if(hd.ncomp > 2) {
                f_ph >> hd.mup[0] >> hd.mup[hd.ncomp-2] >> lnZ0 >> pset;
            } else {
                f_ph >> hd.mup[0] >> hd.mup[hd.ncomp-1] >> lnZ0 >> pset;
            }
            getline(f_ph, dummy);
        } else {
            getline(f_ph, dummy);
            f_ph >> hd.mup[0];
            getline(f_ph, dummy);
        }
        double t, m2min, m2max, m2inc, nmid_val, slope_val;
        while (f_ph >> t >> m2min >> m2max >> m2inc >> nmid_val >> slope_val) {
            t_new.push_back(t);
            mu2min.push_back(m2min);
            mu2max.push_back(m2max);
            mu2_incr.push_back(m2inc);
            hd.nmid.push_back(nmid_val);
            hd.slope.push_back(slope_val);
        }
        ntemp = t_new.size();
        f_ph.close();

    } else if (lpvt) {
        ofstream f_pvt("pvt.dat", ios::app);
        f_pvt << "/* Suffixes and Files = " << suffix;
        for(int ifile=0; ifile<hd.nfiles; ++ifile) f_pvt << hd.filename[ifile] << " ";
        f_pvt << "\n/*    mu1          mu2           <N>           lnZ\n";
        f_pvt.close();

        ifstream f_pvtin("pvt.idat");
        double t;
        while (f_pvtin >> t) {
            t_new.push_back(t);
            vector<double> m_tmp(hd.ncomp);
            for(int icomp=0; icomp<hd.ncomp; ++icomp) f_pvtin >> m_tmp[icomp];
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
            if (mu2_incr[hd.tnum] > 0 && hd.ncomp > 2) hd.mup[hd.ncomp-1] = mu2min[hd.tnum];
            else if (mu2_incr[hd.tnum] < 0 && hd.ncomp > 2) hd.mup[hd.ncomp-1] = mu2max[hd.tnum];
            else if (hd.ncomp == 2) hd.mup[hd.ncomp-1] = mu2min[hd.tnum];

            bool keep_going = true;
            while (keep_going) {
                if (hd.ncomp == 1 && iter == 0) keep_going = true;
                else if (hd.mup[hd.ncomp-1] < mu2max[hd.tnum] + 0.00001 && hd.mup[hd.ncomp-1] > mu2min[hd.tnum] - 0.00001) keep_going = true;
                else keep_going = false;

                if (!keep_going) break;

                if (lfixp) {
                    if (hd.ncomp < 2) { cout << "FIXED PRESSURE REQUIRES NCOMP > 1\n"; return 1; }
                    dmup[hd.ncomp-1] = hd.mu[hd.ncomp-1][0] / 1000.0;
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
                                         (log(nbelow1) - log(nbelow) - log(1.0 - nbelow1) + log(1.0 - nbelow));
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
                            dmup[hd.ncomp-1] /= 2.0;
                            flip2 = 0;
                        }

                        double newsign2;
                        if (sign2 > 0) {
                            hd.mup[hd.ncomp-1] += dmup[hd.ncomp-1];
                            newsign2 = -1;
                        } else {
                            hd.mup[hd.ncomp-1] -= dmup[hd.ncomp-1];
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

                ofstream f_phase("phase.dat", ios::app);
                f_phase << fixed << setprecision(4);
                f_phase << 1.0/hd.betap << " ";
                for(int icomp=0; icomp<hd.ncomp; ++icomp) f_phase << hd.mup[icomp] << " ";
                f_phase << hd.energy_liq << " " << hd.energy_gas << " ";
                for(int icomp=0; icomp<hd.ncomp; ++icomp) f_phase << hd.nliq[icomp] << " ";
                for(int icomp=0; icomp<hd.ncomp; ++icomp) f_phase << hd.ngas[icomp] << " ";
                f_phase << log(hd.Zliq) << " " << log(hd.Zgas) << "\n";
                f_phase.close();

                for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
                    string fname3 = "d" + to_string(icomp+1) + "n" + to_string(hd.tnum+1) + "t" + to_string(num) + "a.dat";
                    ofstream f_d(fname3);
                    f_d << "/* " << fixed << setprecision(3) << 1.0/hd.betap << " ";
                    for(int jcomp=0; jcomp<hd.ncomp; ++jcomp) f_d << fixed << setprecision(2) << hd.mup[jcomp] << " ";
                    f_d << xdim << " " << ydim << " " << zdim << "\n";

                    for (int ipart = hd.m_part[0][icomp]; ipart <= hd.m_part[1][icomp]; ++ipart) {
                        int dens_idx = ipart - hd.minp_global;
                        if(dens_idx >= 0 && dens_idx < hd.dens[icomp].size()) {
                            f_d << setw(6) << ipart << " " << scientific << setprecision(4) << hd.dens[icomp][dens_idx] << "\n";
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
            for (int icomp = 0; icomp < hd.ncomp; ++icomp) avg += hd.avgnum[icomp];

            ofstream f_pvt("pvt.dat", ios::app);
            f_pvt << fixed << setprecision(6) << avg << "     " << log(hd.Z) << "\n";
            f_pvt.close();

            for (int icomp = 0; icomp < hd.ncomp; ++icomp) {
                string fname3 = "d" + to_string(icomp+1) + "n" + to_string(hd.tnum+1) + "t" + to_string(num) + "a.dat";
                ofstream f_d(fname3);
                f_d << "/* " << fixed << setprecision(3) << 1.0/hd.betap << " ";
                for(int jcomp=0; jcomp<hd.ncomp; ++jcomp) f_d << fixed << setprecision(2) << hd.mup[jcomp] << " ";
                f_d << xdim << " " << ydim << " " << zdim << "\n";

                for (int ipart = hd.m_part[0][icomp]; ipart <= hd.m_part[1][icomp]; ++ipart) {
                    int dens_idx = ipart - hd.minp_global;
                    if(dens_idx >= 0 && dens_idx < hd.dens[icomp].size()) {
                        f_d << setw(6) << ipart << " " << scientific << setprecision(4) << hd.dens[icomp][dens_idx] << "\n";
                    }
                }
                f_d.close();
            }
        }
    }

    cout << "Phase coexistence calculation complete\nPROGRAM TERMINATING\n";
    return 0;
}
