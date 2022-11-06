
using PortAudio
#using SampledSignals
using Plots
plotly()
using Base.Threads: @spawn
using Base 
using FFTW
using LinearAlgebra #library used to get the distance between the various points on the map using the norm function
#using PlotlyJS
default(ticks=:native)
default(size=(1080,900))
default(label=""); # Turn off legends\n


function pulse_t(t, pulse_len,pulse_freq,BW)
    A = 0.5 #amplitude of the singal
    rect(t) = (abs.(t) .<= 0.5) * 1.0 #(abs.(t).<=0.5)*1.0; creates a rect function
    sweep_rate = BW / pulse_len
    pulse(t) = A * rect.((t .- pulse_len / 2) / pulse_len) .* cos.(2 * pi * (pulse_freq .* (t .- pulse_len / 2) + 0.5 * sweep_rate * (t .- pulse_len / 2) .^ 2)) #shifting where the tranmitted singal occurs to the centre and the other stretches the signal, the signal now starts at zero and ends at the pulse length.
    #display(plot(pulse(t_Tx)))
    return pulse(t)
end

function filter(pulse_freq, dt, BW, ω)

    A = 1 #amplitude of the singal
    B = 1.2*BW # filter bandwidth in Hz
    # In the sampled frequency domain. add a rect centred on zero to one centred at the next repeat
    # i.e. centred on 0 rad/s an on 2pi/Δt rad/s.
    rect(t) = (abs.(t) .<= 0.5) * 1.0
    H = A * rect((ω .- pulse_freq * (2 * pi)) / (2 * π * B)) .+ rect(((ω .+ pulse_freq * (2 * pi)) .- 2 * π / dt) / (2 * π * B))
    return H
end

function matched_filter_output(v_p) 

    V_T_correlation = FFTW.fft(v_p) #fft of the transmitted signal
    H_MF = conj(V_T_correlation) #creates the matched filter
    return H_MF #cf_RX1#, cf_RX2, cf_RX3
end 


function anayltic_signal(N_samples_R, output_matched_filter_in_frequency)

    v_analytic_frequency =zeros(N_samples_R) #creates an empty array that will hold that data for the analytic signal
    v_analytic_frequency = output_matched_filter_in_frequency #gets the data in the output of the mathced filter and temporarily assigns it to the analytic signal
    v_analytic_frequency[Int(round(N_samples_R/2)):Int(round(N_samples_R))].=0 #sets the seonc half of the array to zero to create the analytic signal
    
    v_analytic_time = (FFTW.ifft(v_analytic_frequency)) #is the analytic signal converted to the time domain
    
    return v_analytic_frequency, v_analytic_time
end


function baseband_output(v_a_t, pulse_freq, t_axis)
    v_baseband_time = v_a_t .* exp.(-im * 2 * pi * pulse_freq .* t_axis) #creates the baseband signal

    return v_baseband_time
end

#displays all the neccessary plots


#function time alsings signals so that they can be subtracted to chekc whether the object is moving or is stationary

function inverse_filter_output(v_p, N_samples_R, dt, pulse_freq, BW, ω)

    B = BW*0.9 # filter bandwidth in Hz

    rect(t) = (abs.(t) .<= 0.5) * 1.0
    H_filtered_so_inverse_can_work = rect((ω .- pulse_freq * (2 * pi)) / (2 * π * B)) .+ rect(((ω .+ pulse_freq * (2 * pi)) .- 2 * π / dt) / (2 * π * B))
    
    V_tx = FFTW.fft(v_p) #this needs to to be the transmitted signal with the correct length

    H_INV=zeros(N_samples_R) #creates an empty array which will store the invrese of the transmitted signal
    H_INV = 1.0./V_tx #inverse of the transmitted signal creating the inverse filter
     #bandlimitng the inverse filter so that no 1/0 points occur
    return H_INV
end

#creates a blackman window for the frequecny points (f) over the bandwidth B
function Blackman(f,B)
    rect(z) = (abs.(z) .<= 0.5) * 1.0
    bm = (0.42 .+ 0.5*cos.((2*pi*f)/B).+0.08*cos.((4*pi*f)/B)).*rect(f/B)
    return bm
end
#function for the the hamming window for the frequecny points (f) over the bandwidth B
function Hamming(f,B)
    rect(z) = (abs.(z) .<= 0.5) * 1.0
    hm = (0.54 .+ 0.46*cos.(2*pi*f/B)).*rect(f/B)
    return hm
end
#function for the the unnamed window for the frequecny points (f) over the bandwidth B
function unnamed_window(f,B)
    rect(z) = (abs.(z) .<= 0.5) * 1.0
    uw = (0.33 .+ 0.66*cos.(pi*f/B).^2).*rect(f/B)
    return uw
end

function time_between_direct_and_reflected(read_current, v_p, t_axis, ω, N_samples_R, dt, c, BW, pulse_freq, f, pulse_len, rx_num, f_values)

    inv_analytic_blackman_time = abs.(FFTW.ifft(window_functions(read_current, v_p, t_axis, ω, N_samples_R, dt, c, BW, pulse_freq, f, pulse_len, rx_num, f_values)))

    first=0
    first_index=0
    second=0
    second_index=0

    for i=2:N_samples_R-1
        current = inv_analytic_blackman_time[i];
        previous= inv_analytic_blackman_time[i-1];
        next= inv_analytic_blackman_time[i+1];
        if  current>previous && current>next && current>first
            second = first;
            second_index = first_index
            first = current;
            first_index= i-1
        elseif  current>previous && current>next && current>second
            second=current
            second_index=i-1
        end
    end

    tau_difference = (second_index-first_index)*dt

    return tau_difference

end
function window_functions(read_current, v_p, t_axis, ω, N_samples_R, dt, c, BW, pulse_freq, f, pulse_len, rx_num, f_values)
    #applies the bandpass filter to the received signal
    read_current_freq = FFTW.fft(read_current).*filter(pulse_freq, dt, BW, ω)

    #applies the inverse filter to the bandpass filtered signal and plots the inverse filtered signal in the time domain
    inverse_read_current_freq = read_current_freq.*inverse_filter_output(v_p, N_samples_R, dt, pulse_freq, BW, ω)

    #applies the mathed filter to the bandpass filtered signal and plots the matched filtered signal
    matched_read_current_freq = read_current_freq.*matched_filter_output(v_p)

    #creates the blackman window function
    window_blackman = Blackman(f .- pulse_freq, BW * 0.9) .+ Blackman((f .+ pulse_freq) .- 1 / dt, BW * 0.9)

    #creates the hamming window function
    window_hamming = Hamming(f .- pulse_freq, BW * 0.9) .+ Hamming((f .+ pulse_freq) .- 1 / dt, BW * 0.9)
    
    #creates the unnamed_window
    window_unamed = unnamed_window(f .- pulse_freq, BW * 0.9) .+ unnamed_window((f .+ pulse_freq) .- 1 / dt, BW * 0.9)
   
    #applies the blackman window to both the inverse and matched filtered signals
    inverse_read_current_freq_window_blackman = inverse_read_current_freq.*window_blackman
    matched_read_current_freq_window_blackman  = matched_read_current_freq.*window_blackman

    #applies the hamming window to both the inverse and matched filtered signals
    inverse_read_current_freq_window_hamming = inverse_read_current_freq.*window_hamming
    matched_read_current_freq_window_hamming  = matched_read_current_freq.*window_hamming

    #applies the unnamed window to both the inverse and matched filtered signals
    inverse_read_current_freq_window_un = inverse_read_current_freq.*window_unamed
    matched_read_current_freq_window_un  = matched_read_current_freq.*window_unamed
   
    #creates the analytic signal for the inverse filtered, blackman windowed function
    inv_analytic_signal_freq_black = inverse_read_current_freq_window_blackman
    inv_analytic_signal_freq_black[Int(round(N_samples_R/2)):N_samples_R].=0

    #creates the analytic signal for the inverse filtered, hamming windowed function
    inv_analytic_signal_freq_ham = inverse_read_current_freq_window_hamming
    inv_analytic_signal_freq_ham[Int(round(N_samples_R/2)):N_samples_R].=0

    #creates the analytic signal for the inverse filtered, unnamed windowed function
    inv_analytic_signal_freq_un = inverse_read_current_freq_window_un
    inv_analytic_signal_freq_un[Int(round(N_samples_R/2)):N_samples_R].=0   
 
    #creates the analytic signal for the matched filtered, blackman windowed function
    match_analytic_signal_freq_black =  matched_read_current_freq_window_blackman
    match_analytic_signal_freq_black[Int(round(N_samples_R/2)):N_samples_R].=0

    #creates the analytic signal for the matched filtered, hamming windowed function
    match_analytic_signal_freq_ham =  matched_read_current_freq_window_hamming
    match_analytic_signal_freq_ham[Int(round(N_samples_R/2)):N_samples_R].=0
    

    #creates the analytic signal for the matched filtered, unnamed windowed function
    match_analytic_signal_freq_un =  matched_read_current_freq_window_un
    match_analytic_signal_freq_un[Int(round(N_samples_R/2)):N_samples_R].=0

    return inv_analytic_signal_freq_black
end

function plot_analytic_signals(t_axis, c,  rx_num, inv_analytic_signal_freq_black, match_analytic_signal_freq_black, inv_analytic_signal_freq_ham, match_analytic_signal_freq_ham, inv_analytic_signal_freq_un, match_analytic_signal_freq_un)
    xlim_min =2
    xlim_max =8

    #creats the plots for both analytic signals with the blakcman window
    plot_inv_ana_black = plot(t_axis*c, abs.(FFTW.ifft(inv_analytic_signal_freq_black)), title="Analytic signal of inverse filtered blackman windowed RX$rx_num", xlims=(xlim_min,xlim_max),  ylabel="Amplitude (V)")
    plot_match_ana_black = plot(t_axis*c, abs.(FFTW.ifft(match_analytic_signal_freq_black)), title="Analytic signal of matched filtered blackman windowed RX$rx_num", xlims=(xlim_min,xlim_max),  xlabel="distance (m)", ylabel="Amplitude (V)")

    #creates the plots for both analytic signals with the hamming window
    plot_inv_ana_ham = plot(t_axis*c, abs.(FFTW.ifft(inv_analytic_signal_freq_ham)), title="Analytic signal of inverse filtered hamming windowed RX$rx_num", xlims=(xlim_min,xlim_max),  ylabel="Amplitude (V)")
    plot_match_ana_ham = plot(t_axis*c, abs.(FFTW.ifft(match_analytic_signal_freq_ham)), title="Analytic signal of matched filtered hamming windowed RX$rx_num", xlims=(xlim_min,xlim_max), ylabel="Amplitude (V)", xlabel="distance (m)")

    #creates the plots for both analytic signals with the hamming window
    plot_inv_ana_un = plot(t_axis*c, abs.(FFTW.ifft(inv_analytic_signal_freq_un)), title="Analytic signal of inverse filtered unnamed windowed RX$rx_num", xlims=(xlim_min,xlim_max),  ylabel="Amplitude (V)")
    plot_match_ana_un = plot(t_axis*c, abs.(FFTW.ifft(match_analytic_signal_freq_un)), title="Analytic signal of matched filtered unnamed windowed RX$rx_num", xlims=(xlim_min,xlim_max), ylabel="Amplitude (V)", xlabel="distance (m)")
    
    #println("works")
    display(plot(plot_inv_ana_black, plot_inv_ana_ham, plot_inv_ana_un, plot_match_ana_black, plot_match_ana_ham, plot_match_ana_un, layout=(3,3)))

end
function locations_of_devices_and_object(Px, Py)
    Tx = [0; 0] #the transimitter position
    P = [Px; Py] #the objects poisiton
    Rx1 = [3; 0] #the first receiver position
    Rx2 = [0; 3] #the second receiver position
    Rx3 = [4; 2] #the third receiver position

    return P, Tx, Rx1, Rx2, Rx3
end

function calculate_distance(P, Tx, Rx1, Rx2, Rx3)
    Tx_P = -Tx + P #vector from Tx to P1

    P_Rx1 = -P + Rx1 #vector from P to Rx1

    Tx_Rx1 = -Tx + Rx1 #vector from Tx to Rx1

    P_Rx2 = -P + Rx2 #vector from P to Rx2

    Tx_Rx2 = -Tx + Rx2 #vector from Tx to Rx2

    P_Rx3 = -P + Rx3 #vector from P to Rx3

    Tx_Rx3 = -Tx + Rx3 #vector from Tx to Rx3

    dist_Tx_P = norm(-Tx + P) #gets euclidean distance from transmitter to P

    dist_P_Rx1 = norm(-P + Rx1) #gets euclidean distance from P to Rx1
    dist_Tx_Rx1 = norm(-Tx + Rx1) # gets euclidean distance from Tx to Rx1

    dist_P_Rx2 = norm(-P + Rx2) #gets euclidean distance from P to Rx2
    dist_Tx_Rx2 = norm(-Tx + Rx2) # gets euclidean distance from Tx to Rx2

    dist_P_Rx3 = norm(-P + Rx3) #gets euclidean distance from P to Rx3
    dist_Tx_Rx3 = norm(-Tx + Rx3) # gets euclidean distance from Tx to Rx3

    return dist_Tx_P, dist_Tx_Rx1, dist_P_Rx1, dist_Tx_Rx2, dist_P_Rx2, dist_Tx_Rx3, dist_P_Rx3
end

function time_delay(dist_Tx_P, dist_Tx_Rx1, dist_P_Rx1, dist_Tx_Rx2, dist_P_Rx2, dist_Tx_Rx3, dist_P_Rx3)

    c = 331 * sqrt(1 + 25 / 273) #speed of light equation that accounts for temperature

    tau_Tx_Rx1 = (dist_Tx_Rx1)/c #time delay from Tx to Rx2 which is known given locations are known

    tau_Tx_P_Rx1 = (dist_P_Rx1 + dist_Tx_P) / c #time delay from the transmitted signal that refelcts of the object and is recieved at Rx1

    #tau_Tx_Rx = dist_Tx_P / c #time delay from the trasnmitter to the reciver

    tau_Tx_Rx2 = (dist_Tx_Rx2)/c #time delay from Tx to Rx2 which is known given locations are known

    tau_Tx_P_Rx2 = (dist_P_Rx2 + dist_Tx_P) / c  #time delay from the transmitted signal that refelcts of the object and is recieved at Rx2

    tau_Tx_Rx3 = (dist_Tx_Rx3)/c #time delay from Tx to Rx3 which is known given locations are known

    tau_Tx_P_Rx3 = (dist_P_Rx3 + dist_Tx_P) / c  #time delay from the transmitted signal that refelcts of the object and is recieved at Rx3

    #tau_Tx_P_Rx_at_TX = (dist_Tx_P+dist_P_Rx_at_Tx)/c; #time delay of signal bouncing off the object and relfecting back to the transmitter where a reciver has been placed.
    return tau_Tx_Rx1, tau_Tx_P_Rx1, tau_Tx_Rx2, tau_Tx_P_Rx2, tau_Tx_Rx3, tau_Tx_P_Rx3
end

function decay_factor(dist_Tx_P, dist_Tx_Rx1, dist_P_Rx1, dist_Tx_Rx2, dist_P_Rx2, dist_Tx_Rx3, dist_P_Rx3)

    A_Tx_Rx1 = 1/(dist_Tx_Rx1^2) #time delay from Tx to Rx2 which is known given locations are known

    A_Tx_P_Rx1 = 1/(dist_P_Rx1^2 * dist_Tx_P^2) #time delay from the transmitted signal that refelcts of the object and is recieved at Rx1

    #tau_Tx_Rx = dist_Tx_P / c #time delay from the trasnmitter to the reciver

    A_Tx_Rx2 = 1/(dist_Tx_Rx2^2) #time delay from Tx to Rx2 which is known given locations are known

    A_Tx_P_Rx2 = 1/(dist_P_Rx2^2 * dist_Tx_P^2)   #time delay from the transmitted signal that refelcts of the object and is recieved at Rx2

    A_Tx_Rx3 = 1/(dist_Tx_Rx3^2) #time delay from Tx to Rx3 which is known given locations are known

    A_Tx_P_Rx3 = 1/(dist_P_Rx3^2 * dist_Tx_P^2)   #time delay from the transmitted signal that refelcts of the object and is recieved at Rx3
 
    return A_Tx_Rx1, A_Tx_P_Rx1, A_Tx_Rx2, A_Tx_P_Rx2, A_Tx_Rx3, A_Tx_P_Rx3
end


function Test_trial_point(TP,TX,RX1,RX2,RX3, tau_Rx1_input, tau_Rx2_input, tau_Rx3_input,ct) #tau_Rxn_input is the tau that i have found
   
    dist_TX_TP = norm(-TX+TP) #temporary distance from tranmitter to trial point
    dist_TP_RX1 = norm(-TP+RX1) #temporary distance from the trial point to Rx1
    dist_TP_RX2 = norm(-TP+RX2) #temporary distance from the trial point to Rx2
    dist_TP_RX3 = norm(-TP+RX3) #temporary distance from the trial point to Rx3

    tau_TX_TP_RX1 = (dist_TX_TP+dist_TP_RX1)/ct #trial point tau value for RX1
    tau_TX_TP_RX2 = (dist_TX_TP+dist_TP_RX2)/ct #trial point tau value for RX2
    tau_TX_TP_RX3 = (dist_TX_TP+dist_TP_RX3)/ct #trial point tau value for RX3

    error_RX1 = (tau_Rx1_input - tau_TX_TP_RX1)^2 #measured tau for Rx1 minus the tau for the trial point from Tx to trial point to RX1 all sqaured
    error_RX2 = (tau_Rx2_input - tau_TX_TP_RX2)^2 #measured tau for Rx2 minus the tau for the trial point from Tx to trial point to RX2 all sqaured
    error_RX3 = (tau_Rx3_input - tau_TX_TP_RX3)^2 #measured tau for Rx3 minus the tau for the trial point from Tx to trial point to RX3 all sqaured

    cost = ((1/2)*(error_RX1+error_RX2+error_RX3))^0.5
end

function find_location(Tx, Rx1, Rx2, Rx3, tau_Tx_P_Rx1, tau_Tx_P_Rx2, tau_Tx_P_Rx3, c, dt)
    #creating a matrix to repreenst the room
    x0 = 0
    y0 = 0
    z0 = 0
    xdim = 4
    ydim = 4
    zdim = 5
    dx = 0.01
    dy = 0.01
    dz = 0.01
    Nx = Int(round(xdim / dx)) + 1
    Ny = Int(round(ydim / dy)) + 1
    Nz = Int(round(zdim / dz)) + 1

    x_axis = x0:dx:(Nx-1)*dx
    y_axis = y0:dy:(Ny-1)*dy

    Store_cost_function = zeros(Nx, Ny)

    #below just runs through the entire 2D map and tests the trial_point in cmparsion to the Tau values found
    for i = 1:Nx
        for j = 1:Ny
            x = x0 + (i - 1) * dx
            y = y0 + (j - 1) * dy
            Trial_point = [x, y]
            Store_cost_function[i, j] = Test_trial_point(Trial_point, Tx, Rx1, Rx2, Rx3, tau_Tx_P_Rx1, tau_Tx_P_Rx2, tau_Tx_P_Rx3, c)

        end
    end
    
    min_point = findmin(Store_cost_function)
    display(surface((Store_cost_function')))

    x_min = (min_point[2][1] - 1) * dx + x0 #subtract one becasue the index starts at 1 and not zero the other math is covnerting from cartesian to how the plot axis are defeined
    y_min = (min_point[2][2] - 1) * dy + y0 #subtract one becasue the index starts at 1 and not zero the other math is covnerting from cartesian to how the plot axis are defeined
    #z_min = (min_point[2][3] - 1) * dz + z0 #subtract one becasue the index starts at 1 and not zero the other math is covnerting from cartesian to how the plot axis are defeined

    return x_min, y_min
end

function main()
    sample_rate = 44100 # sample rate for the signal processing

    dt = 1 / sample_rate #defining the tick rate for the signal time frame

    BW = 2000

    pulse_len = 0.05 #short pulse len

    t_rec = 0.5 #recording time of microphones

    pulse_freq = 13000 #frequency of the pulse

    N_samples_T = Int(round(pulse_len / dt))

    N_samples_R = Int(round(t_rec / dt))

    t_Tx = (0:N_samples_T-1) * dt

    t_axis = 0:dt:(N_samples_R-1)*dt #defines the time axis for plotting

    c = 331 * sqrt(1 + 25 / 273) #speed of light at 25 degrees celcius - assumed tmeperature for expreiment

    lambda = c / pulse_freq #wave length

    Δω = 2 * pi / (N_samples_R * dt) # Sample spacing in freq domain in rad/s
    ω = 0:Δω:(N_samples_R-1)*Δω
    f = ω / (2 * π)

    v_pulse_1 = pulse_t(t_axis, pulse_len, pulse_freq, BW) #transmitted signal


    PSD = 0.01
    N0 = PSD / 2 #for w btn -f0<w<f0
    sigma = sqrt(N0) #variance of the white noise

    #define the location of the object at the start
    Px = 0.8
    Py = 0.5

    println("starting P: ", Px, Py)

    P, Tx, Rx1, Rx2, Rx3 = locations_of_devices_and_object(Px, Py)
    dist_Tx_P, dist_Tx_Rx1, dist_P_Rx1, dist_Tx_Rx2, dist_P_Rx2, dist_Tx_Rx3, dist_P_Rx3 = calculate_distance(P, Tx, Rx1, Rx2, Rx3)
    tau_Tx_Rx1, tau_Tx_P_Rx1, tau_Tx_Rx2, tau_Tx_P_Rx2, tau_Tx_Rx3, tau_Tx_P_Rx3 = time_delay(dist_Tx_P, dist_Tx_Rx1, dist_P_Rx1, dist_Tx_Rx2, dist_P_Rx2, dist_Tx_Rx3, dist_P_Rx3)
    A_Tx_Rx1, A_Tx_P_Rx1, A_Tx_Rx2, A_Tx_P_Rx2, A_Tx_Rx3, A_Tx_P_Rx3 = decay_factor(dist_Tx_P, dist_Tx_Rx1, dist_P_Rx1, dist_Tx_Rx2, dist_P_Rx2, dist_Tx_Rx3, dist_P_Rx3)


    v_Rx1 = (A_Tx_Rx1) .* pulse_t(t_axis .- tau_Tx_Rx1, pulse_len, pulse_freq, BW) .+ (A_Tx_P_Rx1) .* pulse_t(t_axis .- tau_Tx_P_Rx1, pulse_len, pulse_freq, BW) .+ 0.1 * sigma .* randn(N_samples_R) #the recvied pulse, being the transmitted pulse with just a delay from the tranmitter off the object to the Rx1

    v_Rx2 = (A_Tx_Rx2) .* pulse_t(t_axis .- tau_Tx_Rx2, pulse_len, pulse_freq, BW) .+ (A_Tx_P_Rx2) .* pulse_t(t_axis .- tau_Tx_P_Rx2, pulse_len, pulse_freq, BW) .+ 0.1 * sigma .* randn(N_samples_R) #the recvied pulse, being the transmitted pulse with just a delay from the tranmitter off the object to the Rx2

    v_Rx3 = (A_Tx_Rx3) .* pulse_t(t_axis .- tau_Tx_Rx3, pulse_len, pulse_freq, BW) .+ (A_Tx_P_Rx3) .* pulse_t(t_axis .- tau_Tx_P_Rx3, pulse_len, pulse_freq, BW) .+ 0.1 * sigma .* randn(N_samples_R) #the recvied pulse, being the transmitted pulse with just a delay from the tranmitter off the object to the Rx3


    N = length(t_axis)
    if mod(N, 2) == 0    # case N even
        k_values = [0:N/2-1; -N/2:-1]
    else   # case N odd
        k_values = [0:(N-1)/2; -(N - 1)/2:-1]
    end
    df = sample_rate / N
    f_values = k_values * df



    tau_diff_rx1 = time_between_direct_and_reflected(v_Rx1, v_pulse_1, t_axis,  ω, N_samples_R, dt, c, BW, pulse_freq, f, pulse_len, 1, f_values)
    tau_diff_rx2 = time_between_direct_and_reflected(v_Rx2, v_pulse_1, t_axis,  ω, N_samples_R, dt, c, BW, pulse_freq, f, pulse_len, 2, f_values)
    tau_diff_rx3 = time_between_direct_and_reflected(v_Rx3, v_pulse_1, t_axis,  ω, N_samples_R, dt, c, BW, pulse_freq, f, pulse_len, 3, f_values)

    tau_Rx1 = tau_Tx_Rx1 + tau_diff_rx1
    tau_Rx2 = tau_Tx_Rx2 + tau_diff_rx2
    tau_Rx3 = tau_Tx_Rx3 + tau_diff_rx3

    location = find_location(Tx, Rx1, Rx2, Rx3, tau_Rx1, tau_Rx2, tau_Rx3, c, dt)
    println("the object is found to be here: ", location, " the actual place is: ", P)



end

main()