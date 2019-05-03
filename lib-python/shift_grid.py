import numpy as np

NGRID=7

def grid_1d(resol, extent):
    Npix = NGRID*2**resol
    dt = 2*extent/Npix
    grid = np.arange(Npix, dtype=np.float32)*dt + dt/2 - extent
    return grid

def grid_2d(resol, extent):
    x = grid_1d(resol, extent)
    y = grid_1d(resol, extent)
    # convention: x is fast dim, y is slow dim
    grid = np.stack(np.meshgrid(x,y),-1) 
    return grid.reshape(-1,2)

def base_shift_grid(extent):
    return grid_2d(0, extent)

### Neighbor Finding ###

def get_1d_neighbor(mini, curr_res, extent):
    Npix = NGRID*2**(curr_res+1)
    dt = 2*extent/Npix
    ind = np.array([2*mini, 2*mini+1], dtype=np.float32)
    return dt*ind + dt/2 - extent, ind

def get_base_ind(ind):
    xi = ind % NGRID
    yi = ind // NGRID # integer division
    return xi, yi

def get_neighbor(xi, yi, curr_res, extent):
    '''
    Return the 4 nearest neighbors at the next resolution level
    '''
    x_next, xii = get_1d_neighbor(xi, curr_res, extent)
    y_next, yii = get_1d_neighbor(yi, curr_res, extent)
    t_next = np.stack(np.meshgrid(x_next,y_next),-1).reshape(-1,2)
    ind_next = np.stack(np.meshgrid(xii,yii),-1).reshape(-1,2)
    return t_next, ind_next


    

