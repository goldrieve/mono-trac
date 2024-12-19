from geopy.geocoders import Nominatim
from geopy.exc import GeocoderUnavailable
import pandas as pd

# Read the CSV file
locations = pd.read_csv('ml/meta_data.csv')

# Initialize the geolocator
geolocator = Nominatim(user_agent="myGeopyApp")

# Iterate through each entry in the 'Country' column
data = {
    'Country': [],
    'lat': [],
    'lon': []
}

for country in locations['Country']:
    try:
        location = geolocator.geocode(country)
        if location:
            data['Country'].append(country)
            data['lat'].append(location.latitude)
            data['lon'].append(location.longitude)
        else:
            data['Country'].append(country)
            data['lat'].append(None)
            data['lon'].append(None)
    except GeocoderUnavailable:
        data['Country'].append(country)
        data['lat'].append(None)
        data['lon'].append(None)

# Create a pandas DataFrame from the data
df = pd.DataFrame(data)
print(df)