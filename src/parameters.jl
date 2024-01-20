const MANHATTAN_POLY = Tuple{Float32,Float32}[(-74.01369f0,40.69977f0), (-74.00597f0,40.702637f0), (-73.99944f0,40.70641f0), (-73.991714f0,40.708492f0), (-73.9761f0,40.71044f0), (-73.96923f0,40.72931f0), (-73.973526f0,40.736073f0), (-73.9615f0,40.75402f0), (-73.941765f0,40.774693f0), (-73.94348f0,40.78223f0), (-73.938156f0,40.78535f0), (-73.93593f0,40.79029f0), (-73.928894f0,40.79432f0), (-73.92872f0,40.803024f0), (-73.93318f0,40.80744f0), (-73.9349f0,40.833942f0), (-73.92134f0,40.85745f0), (-73.91893f0,40.858356f0), (-73.913956f0,40.863678f0), (-73.909706f0,40.872345f0), (-73.91829f0,40.875168f0), (-73.92648f0,40.879192f0), (-73.93344f0,40.87244f0), (-73.933525f0,40.86793f0), (-73.943436f0,40.853584f0), (-73.947945f0,40.85164f0), (-73.94713f0,40.84414f0), (-73.9552f0,40.828682f0), (-73.96091f0,40.8205f0), (-73.97734f0,40.79864f0), (-73.98957f0,40.78077f0), (-73.996994f0,40.770725f0), (-74.00352f0,40.761368f0), (-74.01064f0,40.75103f0), (-74.01532f0,40.719486f0), (-74.01764f0,40.719063f0), (-74.02047f0,40.704067f0), (-74.01369f0,40.69977f0)];
const CS1 = [(-73.97881, 40.77129), (-73.97107, 40.76796), (-73.97315, 40.76539), (-73.98099, 40.76848), (-73.97881, 40.77129)]
const CS2 = [(-73.97250, 40.78006), (-73.96441, 40.77677), (-73.97010, 40.76916), (-73.97819, 40.77252), (-73.97250, 40.78006)]
const CS3 = [(-73.9690232, 40.7846387), (-73.9710188, 40.7820082), (-73.9705896, 40.7804168), (-73.9639378, 40.7778185), (-73.9624143, 40.7802869), (-73.9645600, 40.7828038), (-73.9690232, 40.7846387)]
const CS4 = [(-73.9638662, 40.7910180),(-73.9679861, 40.7854975),(-73.9633942, 40.7838088),(-73.9614201, 40.7821850),(-73.9570856, 40.7879331),(-73.9638662, 40.7910180)]
const CS5 = [(-73.9577723, 40.7994953),(-73.9632654, 40.7924169),(-73.9553690, 40.7892346),(-73.9502192, 40.7964108),(-73.9577723, 40.7994953)]

############################################################
############################## DATA
############################################################


const TOTX = 9361
const TOTY = 19903
const SEPX = 300  # meters of separation between stations
const SEPY = 300  # meters of separation between stations
const NUMX = ceil(Int, TOTX/SEPX)
const NUMY = ceil(Int, TOTY/SEPY)
const LARGE_NUMBER = 9999999999999
const FULL_MIP_TIME_LIMIT = 10800 # time limit for the full direct model
const MP_TIME_LIMIT = 10800 #10 * 60  # master problem time limit in seconds
const MIP_GAP = 5e-3 # epsilon tolerance for first stage optimality
const GUROBI_ENV = Gurobi.Env()
const NUM_THREADS = 1 # number of threads to be used by GUROBI
const NUM_FOCUS = 0 # default setting for Gurobi parameter NumericFocus
const EPS = 1e-4    # number sufficiently small
const SPEEDS = get_speeds()

const HOR = 45
const OPERATIONAL_HORIZON = 1800
const REFNUM = 2 # one reference stop every REFNUM transit stops
const CAPACITY = 10 # vehicle capacity
const SPEED_FACTOR = 1.
const REF_SPEED_FACTOR = 1.2
const TIME_DISC_MP = 180 # seconds
const TIME_DISC_SP = 30  # seconds
const MAX_DEVIATION_METERS = 750.
const TIME_PERIOD_SECONDS = 900 # time between timeslots/frequencies
const DAY_START = Dates.Time(6) # start of demand horizon at 6 am
const FLEET_SIZE = 10
const THETA = 1.
const MAX_WALKING_METERS = 200.
const MAX_WAITING_SECONDS = 600
const MAX_SCHEDULE_DEV = 600
const WEIGHT_COVERAGE = 100.
const WEIGHT_WALK = 0.01
const WEIGHT_WAIT = 0.01
const WEIGHT_INVEHICLE = 0.01
const WEIGHT_DELAY = 0.01
const WEIGHT_LINE_COST = 0.01
const NORMALIZED = true # noralize values of in-vehicle time and arrival delay with direct taxi trip
const REFSTOP_TIME_GAP = 0 # deviation gap from reference stop time
const EARLY_START_VEHICLES = 1800 # how earlier in seconds than demand horizon (DAY_START) should first vehicle start operating (operating horizon)
const WALKSPEED_M_PER_SEC = 1.4

# exit points toward Laguardia: 1. Kennedy bridge (7 min),  2. Queensboro bridge (18 min), 3. Williamsburg (14 min), 4. Midtown tunnel (20 min)
# const EXIT_POINTS_OSM = Int[point_to_nodes(LLA(40.80267437304975, -73.93345544007066), m), point_to_nodes(LLA(40.76090457998287, -73.96415881726622), m), point_to_nodes(LLA(40.7182818758882, -73.9873787003149), m), point_to_nodes(LLA(40.74703611573972, -73.9770968903696), m)]
const EXIT_POINTS_OSM = Int[42423203, 1775176473, 486867534, 42449932]
const EXIT_TIMES_SECONDS = Float64[7*60, 18*60, 14*60, 20*60]
const M_PER_MILE = 1609.34
const EXIT_DIST_METERS = Float64[5.3*M_PER_MILE, 7.4*M_PER_MILE, 8.1*M_PER_MILE, 7.2*M_PER_MILE]

# additional MiND-DAR parameters
const REF_SPEED_FACTOR_DAR = 1.1
const TIME_DISC_MP_DAR = 300 # seconds