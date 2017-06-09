import netCDF4 as nc
import numpy as np

class StatsFile:
    def __init__(self, file):
        self.file = file

        return


    def mean(self, var_name, group, tmin=None, tmax=None, zmin=None, zmax=None):


        if group == 'timeseries':
            return self.__mean_ts(var_name, tmin, tmax)
        elif group == 'profiles':
            return self.__mean_ps(var_name, tmin, tmax, zmin, zmax)



        return

    def data(self, var_name, group, tmin=None, tmax=None, zmin=None, zmax=None):



        return


    def __mean_ts(self, var_name, tmin=None, tmax=None):

        rt_grp = nc.Dataset(self.file, 'r')
        ts = rt_grp['timeseries']

        time = ts['t'][:]


        if tmin is None:
            t_min_idx = 0
        else:
            t_min_idx = np.where(np.abs(time - tmin) == np.min(np.abs(time - tmin)))[0][0]


        if tmax is None:
            t_max_idx = -1
        else:
            t_max_idx = np.where(np.abs(time - tmax) == np.min(np.abs(time - tmax)))[0][0]


        mean = np.mean(ts[var_name][t_min_idx:t_max_idx])

        rt_grp.close()

        return mean


    def __mean_ps(self, var_name, tmin=None, tmax=None, zmin=None, zmax=None):


        rt_grp = nc.Dataset(self.file, 'r')
        ps = rt_grp['profiles']
        ref = rt_grp['reference']


        time = ps['t'][:]
        z = ref['zp_half'][:]

        if tmin is None:
            t_min_idx = 0
        else:
            t_min_idx = np.where(np.abs(time - tmin) == np.min(np.abs(time - tmin)))[0][0]


        if tmax is None:
            t_max_idx = -1
        else:
            t_max_idx = np.where(np.abs(time - tmax) == np.min(np.abs(time - tmax)))[0][0]




        if zmin is None:
            z_min_idx = 0
        else:
            z_min_idx = np.where(np.abs(z - zmin) == np.min(np.abs(z - zmin)))[0][0]

        if zmax is None:
            z_max_idx = -1
        else:
            z_max_idx = np.where(np.abs(z - zmax) == np.min(np.abs(z - zmax)))[0][0]



        mean = np.mean(ps[var_name][t_min_idx:t_max_idx,z_min_idx:z_max_idx],axis=0)

        rt_grp.close()


        print np.shape(mean)
        return mean, z[z_min_idx:z_max_idx]



def main():


    StatsFile1 = StatsFile('/Users/kpressel/Desktop/RicoData/SB/stats/Stats.Rico.nc')
    StatsFile2 = StatsFile('/Users/kpressel/Desktop/RicoData/CoSal/stats/Stats.Rico.nc')
    StatsFile3 = StatsFile('/Users/kpressel/Desktop/RicoData/Co70/stats/Stats.Rico.nc')

    print  StatsFile2.mean('lwp', 'timeseries', tmin=86400-1*3600, tmax=86400)

    import pylab as plt
    plt.figure()
    mean, z = StatsFile1.mean('qr_mean', 'profiles', tmin=86400-4*3600, tmax=86400)
    plt.plot(mean*1000000.0, z)
    mean, z = StatsFile2.mean('qrain_mean', 'profiles', tmin=86400-4*3600, tmax=86400)
    plt.plot(mean*1000000.0, z)
    mean, z = StatsFile3.mean('qrain_mean', 'profiles', tmin=86400-4*3600, tmax=86400)
    plt.plot(mean*1000000.0, z)
    plt.grid()
    plt.show()


    return



if __name__ == '__main__':
    main()