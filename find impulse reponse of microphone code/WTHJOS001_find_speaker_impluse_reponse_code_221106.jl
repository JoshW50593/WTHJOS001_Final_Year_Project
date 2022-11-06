cd("D:\\UCT 2022\\term 2\\4022\\testing_mic")

using PortAudio
using SampledSignals
using Plots
#plotly()
gr()
default(size=(1080,900))

using Base 
using Base.Threads: @spawn
using WAV
using FFTW
using LinearAlgebra #library used to get the distance between the various points on the map using the norm function
#using PlotlyJS

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
    B = 1.1*BW # filter bandwidth in Hz
    # In the sampled frequency domain. add a rect centred on zero to one centred at the next repeat
    # i.e. centred on 0 rad/s an on 2pi/Δt rad/s.
    rect(t) = (abs.(t) .<= 0.5) * 1.0
    H = A * rect((ω .- pulse_freq * (2 * pi)) / (2 * π * B)) .+ A* rect((ω .+ pulse_freq * (2 * pi).- 2 * pi / dt) / (2 * π * B))
    return H
end
#opens the stream to a usb port on the compture to a micrphone
function read_streams_2(sample_rate, N_samples_R)
    stream_2 =PortAudio.PortAudioStream("Microphone (4- USB PnP Sound De",1,0, samplerate=sample_rate) #change the micrphone name for the different micrphones used
    v_Rx_2 = read(stream_2,N_samples_R) #gets the data receive by the microphone
    ans_2=v_Rx_2
    close(stream_2) #closes the read microphone stream
    return ans_2
end

#generates a transmitted pulse with n pulses
function v_tx(n, t_axis, pulse_len, pulse_freq, BW, N_samples_R)
    pulse_shift_factor =50
    v_pulse=zeros(N_samples_R)
    for i in 1:Int(n)
        v_pulse+=pulse_t(t_axis.-(pulse_len*pulse_shift_factor*i), pulse_len, pulse_freq, BW)
    end
    return v_pulse
end

function matched_filter_output(v_p) # V_rx2_filtered, V_rx3_filtered,


    V_T_correlation = FFTW.fft(v_p) #fft of the transmitted signal

    H_MF = conj(V_T_correlation) #creates the matched filter

    return H_MF #cf_RX1#, cf_RX2, cf_RX3
end 



function main()
    sleep(1)
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

    n = 7 #the number of pulses in the transmitted signal

    v_Tx = v_tx(n, t_axis, pulse_len, pulse_freq, BW, N_samples_R)
    global v_pulse = pulse_t(t_axis, pulse_len, pulse_freq, BW) #propagated singal that is not yet truncated
 

    N = length(t_axis)


    if mod(N, 2) == 0    # case N even
        k_values = [0:N/2-1; -N/2:-1]
    else   # case N odd
        k_values = [0:(N-1)/2; -(N - 1)/2:-1]
    end
    df = sample_rate / N

    f_values = k_values * df
  
    #writes to stream the output stream and reads from a given input stream
    out_stream = PortAudio.PortAudioStream(0, 2) #opens teh output stream to the speaker
    global v_stream_2 = @spawn read_streams_2(sample_rate, N_samples_R) #reads from stream two
    write(out_stream, v_pulse) #writes the pulse to the output stream
    close(out_stream) #closes the output stream
    v_Rx2 = Base.fetch(v_stream_2) #gets the mircphone data and stores it in v_Rx2

    #checks for clipping of the signal
    if maximum(Base.fetch(v_stream_2))>0.98 || minimum(Base.fetch(v_stream_2))<-0.98
        println("clipping")
    end

    #second read


    tf = (FFTW.fft(v_Rx2).*filter(pulse_freq, dt, BW, ω))./FFTW.fft(v_pulse)
    #used to write the impulse response for one spekaer-micrphone to a wav file
    wavwrite(real(FFTW.ifft(tf)), "impulse_rx_time_mic_with_tape_over_used_freq_221102.wav", Fs=44100)
    #when a different microphone is used ensure to change the wav file name to somethin suitable or just ensure it is not the same as teh previous one otherwise it will be overwritten
 

end

main()