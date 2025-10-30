using Lambert
using Plots
using AstroCoords
using SciMLBase

# Sun's gravitational parameter
μ_sun = 1.32712440018e11  # [km³/s²]

# Earth orbital elements (osculating at J2000 epoch)
# Standard reference values from JPL HORIZONS
earth_a = 149.598023e6      # Semi-major axis [km]
earth_e = 0.0167086         # Eccentricity
earth_i = 0.00005 * π/180   # Inclination [rad] (essentially 0° - ecliptic reference)
earth_Ω = 0.0               # RAAN [rad]
earth_ω = 102.94719 * π/180 # Argument of periapsis [rad]
earth_ν0 = 0.0              # Initial true anomaly [rad] (start at perihelion)

# Create Earth Keplerian state
earth_kep0 = Keplerian(earth_a, earth_e, earth_i, earth_Ω, earth_ω, earth_ν0)

# Mars orbital elements (osculating at J2000 epoch)
# Standard reference values from JPL HORIZONS
mars_a = 227.923722e6       # Semi-major axis [km]
mars_e = 0.0934             # Eccentricity
mars_i = 1.850 * π/180      # Inclination [rad]
mars_Ω = 49.558 * π/180     # RAAN [rad]
mars_ω = 286.502 * π/180    # Argument of periapsis [rad]
mars_ν0 = 19.412 * π/180    # Initial true anomaly [rad] (phase offset from Earth)

# Create Mars Keplerian state
mars_kep0 = Keplerian(mars_a, mars_e, mars_i, mars_Ω, mars_ω, mars_ν0)

# Calculate orbital periods using AstroCoords
earth_period = orbitalPeriod(earth_kep0, μ_sun)
mars_period = orbitalPeriod(mars_kep0, μ_sun)

println("Earth orbital period: $(round(earth_period / 86400, digits=2)) days")
println("Mars orbital period: $(round(mars_period / 86400, digits=2)) days")

# Function to propagate Keplerian orbit using mean anomaly
function propagate_orbit(kep0::Keplerian, t::Number, μ::Number)
    # Get mean motion from AstroCoords
    n = meanMotion(kep0, μ)

    # Update mean longitude
    M_new = kep0.M + n * t

    # Create new Keplerian with updated mean longitude
    kep_new = Keplerian(kep0.a, kep0.e, kep0.i, kep0.Ω, kep0.ω, M_new)

    # Convert to Cartesian to get position and velocity
    cart = Cartesian(kep_new, μ)

    return AstroCoords.params(cart)
end

# Create state propagation functions (return [r; v])
state_earth(t::Number) = propagate_orbit(earth_kep0, t, μ_sun)
state_mars(t::Number) = propagate_orbit(mars_kep0, t, μ_sun)

# Define departure and arrival windows
days_to_seconds = 86400.0
departure_window = range(0, 600 * days_to_seconds, length = 100)  # 0 to 200 days (in seconds)
arrival_window = range(300 * days_to_seconds, 1200 * days_to_seconds, length = 100)  # 150 to 350 days (in seconds)

println("\nCreating Earth-Mars transfer porkchop plot...")
println(
    "This will evaluate $(length(departure_window) * length(arrival_window)) Lambert problems...",
)
println("This may take a minute or two...")

p = porkchop_plot(
    μ_sun,
    state_earth,
    state_mars,
    departure_window,
    arrival_window;
    time_scale = days_to_seconds,  # Scale axes to show days instead of seconds
    max_deltav = 30.0,     # Filter trajectories requiring > 30 km/s total ΔV
    title = "Earth to Mars Transfer - Porkchop Plot",
    xlabel = "Departure Time [days]",
    ylabel = "Arrival Time [days]",
    tof_spacing = 200 * days_to_seconds,
    plot_quantity = :total_excess_velocity,
)

println("✓ Porkchop plot created successfully!")
println("  Plot type: ", typeof(p))
println("  Number of series: ", length(p.series_list))

println("\nSaving plot to 'earth_mars_porkchop_total_excess_velocity.png'...")
savefig(p, "earth_mars_porkchop_total_excess_velocity.png")
println("✓ Plot saved!")

println("\n✓ Earth-Mars transfer porkchop plot complete!")
println("Look for the characteristic 'porkchop' shape showing optimal transfer windows.")
println("The minimum ΔV regions correspond to Hohmann-like transfers.")

p2 = porkchop_plot(
    μ_sun,
    state_earth,
    state_mars,
    departure_window,
    arrival_window;
    time_scale = days_to_seconds,  # Scale axes to show days instead of seconds
    max_deltav = 30.0,     # Filter trajectories requiring > 30 km/s total ΔV
    title = "Earth to Mars Transfer - Porkchop Plot",
    xlabel = "Departure Time [days]",
    ylabel = "Arrival Time [days]",
    tof_spacing = 200 * days_to_seconds,
    plot_quantity = :excess_velocity_departure,
)

println("✓ Porkchop plot created successfully!")
println("  Plot type: ", typeof(p2))
println("  Number of series: ", length(p2.series_list))

println("\nSaving plot to 'earth_mars_porkchop_excess_velocity_departure.png'...")
savefig(p2, "earth_mars_porkchop_excess_velocity_departure.png")
println("✓ Plot saved!")

println("\n✓ Earth-Mars transfer porkchop plot complete!")
println("Look for the characteristic 'porkchop' shape showing optimal transfer windows.")
println("The minimum ΔV regions correspond to Hohmann-like transfers.")

p3 = porkchop_plot(
    μ_sun,
    state_earth,
    state_mars,
    departure_window,
    arrival_window;
    time_scale = days_to_seconds,  # Scale axes to show days instead of seconds
    max_deltav = 30.0,     # Filter trajectories requiring > 30 km/s total ΔV
    title = "Earth to Mars Transfer - Porkchop Plot",
    xlabel = "Departure Time [days]",
    ylabel = "Arrival Time [days]",
    tof_spacing = 200 * days_to_seconds,
    plot_quantity = :excess_velocity_arrival,
)

println("✓ Porkchop plot created successfully!")
println("  Plot type: ", typeof(p3))
println("  Number of series: ", length(p3.series_list))

println("\nSaving plot to 'earth_mars_porkchop_excess_velocity_arrival.png'...")
savefig(p3, "earth_mars_porkchop_excess_velocity_arrival.png")
println("✓ Plot saved!")

println("\n✓ Earth-Mars transfer porkchop plot complete!")
println("Look for the characteristic 'porkchop' shape showing optimal transfer windows.")
println("The minimum ΔV regions correspond to Hohmann-like transfers.")


p4 = porkchop_plot(
    μ_sun,
    state_earth,
    state_mars,
    departure_window,
    arrival_window;
    time_scale = days_to_seconds,  # Scale axes to show days instead of seconds
    max_deltav = 30.0,     # Filter trajectories requiring > 30 km/s total ΔV
    title = "Earth to Mars Transfer - Porkchop Plot",
    xlabel = "Departure Time [days]",
    ylabel = "Arrival Time [days]",
    tof_spacing = 200 * days_to_seconds,
    plot_quantity = [
        :total_excess_velocity,
        :excess_velocity_departure,
        :excess_velocity_arrival,
    ],
)

println("✓ Porkchop plot created successfully!")
println("  Plot type: ", typeof(p4))
println("  Number of series: ", length(p4.series_list))

println("\nSaving plot to 'earth_mars_porkchop_all.png'...")
savefig(p4, "earth_mars_porkchop_all.png")
println("✓ Plot saved!")
