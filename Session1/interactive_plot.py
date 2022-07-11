from matplotlib.widgets import Slider, Button
import numpy as np
import matplotlib.pyplot as plt
def calcW_n(k:float,m:float):
    """
    Implement the expression of W_n

    Arguments:
    k, m -- parameters of the SDOF system

    Returns:
    w_n -- the value of w_n, computed using the formula seen in the course
    """   

    ### BEGIN SOLUTION
    w_n = np.sqrt(k/m)
    ### END SOLUTION

    return w_n

def calcXi(c:float,m:float,k:float):
    """
    Implement the expression of W_n
    
    Arguments:
    k,m -- parameter of the SDOF system
    
    Returns:
    w_n -- the value of w_n, computed using the formula seen in the course
    """   

    ### BEGIN SOLUTION
    Xi= c/(2*np.sqrt(m*k))
    ### END SOLUTION
    return Xi

def impulse_response(m:float,c:float,k:float):
    """
    Implement the SDOF impulse response 
    
    Arguments:
    m, c and k -- Parameter of the system 
                -- m is in kg
                -- c is in N.s.m^-1
                -- k is in N.m^-1
    
    Returns:
    h -- the value of function h, computed using 
        expression of the impulse response
    -----------------------------
    tips :  1.call the function defined above to compute Xi and Wn and 
              save them in a local variable
            2.define a variable wd first
            3.use np.exp and np.sin
    -----------------------------
    """
    # as we are computing a numerical solution we need to define 
    # a time vector
    t =np.arange(start=0,stop=100,step=0.1)
    
    ### BEGIN SOLUTION
    Xi=calcXi(c=c,m=m,k=k)
    w_n=calcW_n(k=k,m=m)
    
    w_d= w_n*np.sqrt(1-Xi**2)
    h=1/(m*w_d)*np.exp(-Xi*w_d*t)*np.sin(w_d*t)
    ### END SOLUTION
    return h ,t
    
def SDOF_plot():
    # Define initial parameters
    init_m = 1
    init_k = 16
    init_c = 0.1

    # Create the figure and the line that we will manipulate
    fig, ax = plt.subplots()

    h,t=impulse_response(m=init_m,k=init_k,c=init_c)
    line, = plt.plot(t,h, lw=2)
    ax.set_xlabel('Time (s)')

    # adjust the main plot to make room for the sliders
    plt.subplots_adjust(left=0.25, bottom=0.4)

    wn=calcW_n(k=init_k,m=init_m)
    Xi=calcXi(c=init_c,k=init_k,m=init_m)
    wd= wn*np.sqrt(1-Xi**2)

    text_wn = plt.text(-43, 0.25, f'$\omega_n$ = {wn:.2f}',
            style ='italic',
            fontsize = 10,
            bbox ={'facecolor':'green',
                'alpha':0.6, 'pad':2})

    text_Xi = plt.text(-43, 0.2, f'$ \\xi $ = {Xi*100:.2f}%',
            style ='italic',
            fontsize = 10,
            bbox ={'facecolor':'green',
                'alpha':0.6, 'pad':2})

    text_wd = plt.text(-43, 0.15, f'$ \omega_d $ = {wd:.2f}',
            style ='italic',
            fontsize = 10,
            bbox ={'facecolor':'green',
                'alpha':0.6, 'pad':2})

    axmass = plt.axes([0.25, 0.1, 0.65, 0.03])
    m_slider = Slider(
        ax=axmass,
        label='Mass (kg)',
        valmin=0.1,
        valmax=30,
        valinit=init_m,
    )
    axstifness=plt.axes([0.25, 0.15, 0.65, 0.03])
    k_slider = Slider(
        ax=axstifness,
        label="Stifness (N/m)",
        valmin=0.1,
        valmax=20,
        valinit=init_k,
    )
    axc = plt.axes([0.25, 0.2, 0.65, 0.03])
    c_slider = Slider(
        ax=axc,
        label="friction (N/(m*s))",
        valmin=0.01,
        valmax=10,
        valinit=init_c,
    )



    # The function to be called anytime a slider's value changes
    def update():
        h,t=impulse_response(m=m_slider.val,k=k_slider.val,c=c_slider.val)
        line.set_ydata(h)
        Xi=calcXi(c=c_slider.val,k=k_slider.val,m=m_slider.val)
        wn=calcW_n(k=k_slider.val,m=m_slider.val)
        wd= wn*np.sqrt(1-Xi**2)

        text_wn.set_text(f'$\omega_n$ = {wn:.2f}')
        text_Xi.set_text(f'$ \\xi $ = {Xi*100:.2f}%')
        text_wd.set_text(f'$ \omega_d $ = {wd:.2f}')
        
        fig.canvas.draw_idle()

    # register the update function with each slider
    m_slider.on_changed(update)
    k_slider.on_changed(update)
    c_slider.on_changed(update)


    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')


    def reset(event):
        m_slider.reset()
        k_slider.reset()
    button.on_clicked(reset)

    return update