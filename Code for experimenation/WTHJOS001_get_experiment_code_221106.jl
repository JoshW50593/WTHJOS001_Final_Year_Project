
using PortAudio
#using SampledSignals
using Plots
pyplot()
using Base 
using Base.Threads: @spawn
using WAV
using FFTW
using LinearAlgebra #library used to get the distance between the various points on the map using the norm function
#using PlotlyJS

function pulse_t(t, pulse_len,pulse_freq,BW)
    A = 0.5 #amplitude of the singal
    rect(t) = (abs.(t) .<= 0.5) * 1.0 #(abs.(t).<=0.5)*1.0; creates a rect function
    #BW = 5000
    sweep_rate = BW / pulse_len
    pulse(t) = A * rect.((t .- pulse_len / 2) / pulse_len) .* cos.(2 * pi * (pulse_freq .* (t .- pulse_len / 2) + 0.5 * sweep_rate * (t .- pulse_len / 2) .^ 2)) #shifting where the tranmitted singal occurs to the centre and the other stretches the signal, the signal now starts at zero and ends at the pulse length.
    #display(plot(pulse(t_Tx)))
    return pulse(t)
end

function read_streams_1(sample_rate, N_samples_R) #tape
    stream = PortAudio.PortAudioStream("Microphone (2- USB PnP Sound De",1,0, samplerate=sample_rate)
    v_Rx_1 = read(stream,N_samples_R)
    ans_1=v_Rx_1
    return ans_1
end

function read_streams_2(sample_rate, N_samples_R) #nothing
    stream_2 =PortAudio.PortAudioStream("Microphone (4- USB PnP Sound De",1,0, samplerate=sample_rate) 
    v_Rx_2 = read(stream_2,N_samples_R)
    ans_2=v_Rx_2
    return ans_2
end

function read_streams_3(sample_rate, N_samples_R) #prestick
    stream_3 = PortAudio.PortAudioStream("Microphone (USB PnP Sound Devic",1,0,samplerate=sample_rate)
    v_Rx_3 = read(stream_3, N_samples_R)
    ans_3 = v_Rx_3
    return ans_3
end

function write_to_stream(v_p)
    out_stream = PortAudio.PortAudioStream(0,2)
    @spawn write(out_stream,v_p)
    ans = v_p
    return ans
end

function filter(pulse_freq, dt, BW, ω)

    A = 1 #amplitude of the singal
    B = 1.1*BW # filter bandwidth in Hz
    # In the sampled frequency domain. add a rect centred on zero to one centred at the next repeat
    # i.e. centred on 0 rad/s an on 2pi/Δt rad/s.
    rect(t) = (abs.(t) .<= 0.5) * 1.0
    H = A * rect((ω .- pulse_freq * (2 * pi)) / (2 * π * B)) .+ rect(((ω .+ pulse_freq * (2 * pi)) .- 2 * π / dt) / (2 * π * B))
    return H
end

function matched_filter_output(v_p) # V_rx2_filtered, V_rx3_filtered,

    V_T_correlation = FFTW.fft(v_p) #fft of the transmitted signal

    H_MF = 1.0./ V_T_correlation #creates the matched filter
    return H_MF #cf_RX1#, cf_RX2, cf_RX3
end 



function Blackman(f, B)
    rect(z) = (abs.(z) .<= 0.5) * 1.0
    bm = (0.42 .+ 0.5 * cos.((2 * pi * f) / B) .+ 0.08 * cos.((4 * pi * f) / B)) .* rect(f / B)
    return bm
end
#displays all the neccessary plot

#function time alsings signals so that they can be subtracted to chekc whether the object is moving or is stationary


function locations_of_devices_and_object(Px,Py,Pz,distance_moved_x,distance_moved_y,distance_moved_z)
    Tx = [0; 0] #the transimitter position
    #Rx_at_Tx=[.1;0;0] #for the intial simuation purposes we will assume that this is bascially on top of the transmitter
    Rx2 = [0.5; 1.5] #the recieve position
    P = [Px+distance_moved_x; Py+distance_moved_y] #the objects poisiton
    Rx3 = [1.5; 0.5]
    Rx1 = [0; 0.5]
    P_initialized =[Px;Py]
    #P_moved_only = [Px+distance_moved_x;Py+distance_moved_y;Pz+distance_moved_z]
    #println(P," location of object that has moved")
    return P, Tx, Rx1, Rx2, Rx3

end


function Blackman(f, B)
    rect(z) = (abs.(z) .<= 0.5) * 1.0
    bm = (0.42 .+ 0.5 * cos.((2 * pi * f) / B) .+ 0.08 * cos.((4 * pi * f) / B)) .* rect(f / B)
    return bm
end

function Hamming(f, B)
    rect(z) = (abs.(z) .<= 0.5) * 1.0
    bm = (0.54.+ 0.46 * cos.((2 * pi * f) / B)) .* rect(f / B)
    return bm
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

    pulse_shift_factor =4
    v_pulse_1 = pulse_t(t_axis, pulse_len, pulse_freq, BW) #propagated singal that is not yet truncated
    v_pulse_2 = pulse_t(t_axis.-(pulse_len*pulse_shift_factor), pulse_len, pulse_freq, BW)
    v_pulse_3 = pulse_t(t_axis.-(pulse_len*2*pulse_shift_factor), pulse_len, pulse_freq, BW)
    v_pulse_4 = pulse_t(t_axis.-(pulse_len*3*pulse_shift_factor), pulse_len, pulse_freq, BW)
    v_pulse_5 = pulse_t(t_axis.-(pulse_len*4*pulse_shift_factor), pulse_len, pulse_freq, BW)

    v_pulse = v_pulse_1#.+v_pulse_2.+v_pulse_3.+v_pulse_4.+v_pulse_5
    #the PSD for white nosie can be classifed by being equal N0/2 for frewuency within 2B cntered at f=0


    out_stream = PortAudio.PortAudioStream(0, 2)

    v_stream_2 = @spawn read_streams_2(sample_rate, N_samples_R) #nothing
    v_stream_3 = @spawn read_streams_3(sample_rate, N_samples_R) #prestick
    v_stream_1 = @spawn read_streams_1(sample_rate, N_samples_R)# tape
    write(out_stream, v_pulse)

    close(out_stream)

    ref_impulse_prestick = wavread("impulse_rx_time_mic_with_prestick_over_used_freq_221102.wav")
    ref_impulse_nothing = wavread("impulse_rx_time_mic_with_nothing_over_used_freq_221102.wav")
    ref_impulse_tape = wavread("impulse_rx_time_mic_with_tape_over_used_freq_221102.wav")

    global tf_prestick = (FFTW.fft(ref_impulse_prestick[1]))
    global tf_nothing = (FFTW.fft(ref_impulse_nothing[1]))
    global tf_tape = (FFTW.fft(ref_impulse_tape[1]))


    sleep(2)

    for i = 1:10

        out_stream = PortAudio.PortAudioStream(0, 2)
        global v_Rx2 = @spawn read_streams_2(sample_rate, N_samples_R) #nothing
        global v_Rx3 = @spawn read_streams_3(sample_rate, N_samples_R) #prestick
        global v_Rx1 = @spawn read_streams_1(sample_rate, N_samples_R) #tape
        global v_p = @spawn write_to_stream(v_pulse)
        
        read_current_Rx1 = FFTW.fft(Base.fetch(v_Rx1)).*filter(pulse_freq, dt, BW, ω)
        read_current_Rx2 = FFTW.fft(Base.fetch(v_Rx2)).*filter(pulse_freq, dt, BW, ω)
        read_current_Rx3 = FFTW.fft(Base.fetch(v_Rx3)).*filter(pulse_freq, dt, BW, ω)

        wavwrite(Base.fetch(v_Rx1), "v_rx1_sampled_signal_mulitstatic_1_221102_$i.wav", Fs=44100)
        wavwrite(Base.fetch(v_Rx2), "v_rx2_sampled_signal_multistatic_1_221102_$i.wav", Fs=44100)
        wavwrite(Base.fetch(v_Rx3), "v_rx3_sampled_signal_multistatic_1_221022_$i.wav", Fs=44100)

        vrx_tf_tape = ((read_current_Rx1)./(tf_tape)).*filter(pulse_freq, dt, BW, ω)
        vrx_tf_nothing = ((read_current_Rx2)./(tf_nothing)).*filter(pulse_freq, dt, BW, ω)
        vrx_tf_prestick = ((read_current_Rx3)./(tf_prestick)).*filter(pulse_freq, dt, BW, ω)


        plot_Rx1_current = plot(t_axis*c, real(FFTW.ifft(read_current_Rx1)), title="rx1 filtered")
    

        inverse_tape = vrx_tf_tape.*matched_filter_output(v_pulse)
        inverse_prestick = vrx_tf_prestick.*matched_filter_output(v_pulse)
        inverse_nothing = vrx_tf_nothing.*matched_filter_output(v_pulse)

        inverse_tape[Int(round(N_samples_R/2)):Int(round(N_samples_R))].=0
        inverse_prestick[Int(round(N_samples_R/2)):Int(round(N_samples_R))].=0
        inverse_nothing[Int(round(N_samples_R/2)):Int(round(N_samples_R))].=0

        window = Blackman(f .- pulse_freq, BW * 0.9) .+ Blackman((f .+ pulse_freq) .- 1 / dt, BW * 0.9)

        blackman_inverse_tape = inverse_tape.*window
        blackman_inverse_prestick = inverse_prestick.*window
        blackman_inverse_nothing = inverse_nothing.*window
        #divided_rx1 = ((read_current_Rx1) ./ tf).*filter(pulse_freq, dt, BW, ω)
        plot_divided_tape = plot(t_axis*c, abs.(ifft(blackman_inverse_tape)), title="Analytic signal Rx1", xlims=(15, 70))
        plot_divided_prestick = plot(t_axis*c, abs.(ifft(blackman_inverse_prestick)), title="Analytic signal Rx2", xlims=(15,70))
        plot_divided_nothing = plot(t_axis*c, abs.(ifft(blackman_inverse_nothing)), title="Anayltic signal Rx3", xlims=(15,70))
        #matched = (divided_rx1).*matched_filter_output(v_pulse)

        inv_analytic_blackman_time_tape = abs.(FFTW.ifft(blackman_inverse_tape ))
        inv_analytic_blackman_time_prestick = abs.(FFTW.ifft(blackman_inverse_prestick))
        inv_analytic_blackman_time_nothing = abs.(FFTW.ifft(blackman_inverse_nothing))

        wavwrite(inv_analytic_blackman_time_tape, "output_of_inverse_filter_with_blakcman_filter_for_mic_with_tape_in_time_domain_$i.wav", Fs=44100)
        wavwrite(inv_analytic_blackman_time_prestick, "output_of_inverse_filter_with_blakcman_filter_for_mic_with_prestick_in_time_domain_$i.wav", Fs=44100)
        wavwrite(inv_analytic_blackman_time_nothing, "output_of_inverse_filter_with_blakcman_filter_for_mic_with_nothing_in_time_domain_$i.wav", Fs=44100)
        #matched_plot = plot(t_axis*c, abs.(FFTW.ifft(matched)), title="matched")
        pyplot()
        display(plot(plot_divided_tape, plot_divided_prestick, plot_divided_nothing, layout=(3,1)))
        sleep(30)
    end

end

main()