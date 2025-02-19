from geopy.geocoders import Nominatim
import pandas as pd
from geopy.extra.rate_limiter import RateLimiter

# Read the CSV file
locations = pd.read_csv('ml/meta_data.csv')

# Initialize the geolocator
geolocator = Nominatim(user_agent="application")


reverse = RateLimiter(geolocator.reverse, min_delay_seconds=1)

location = reverse((50.6539239, -120.3385242), language='en', exactly_one=True)

print(location.raw)