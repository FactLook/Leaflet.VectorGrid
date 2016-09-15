#!/usr/bin/env bash

# Helper script to create browser-loadable versions of geojson-vt, pbf and vector-tile
cd ..

cd node_modules/geojson-vt/
npm install
npm run build-dev
cd -

cd node_modules/pbf/
npm install
npm run build-dev
cd -

cd node_modules/vector-tile/
npm install
mkdir -p dist && ../../node_modules/pbf/node_modules/.bin/browserify index.js -d -s vector-tile > dist/vector-tile-dev.js
cd -

cd node_modules/leaflet/
npm install
jake build
cd -

mkdir -p vendor
cp node_modules/geojson-vt/geojson-vt-dev.js vendor/
cp node_modules/pbf/dist/pbf-dev.js vendor/
cp node_modules/vector-tile/dist/vector-tile-dev.js vendor/
cp node_modules/leaflet/dist/leaflet* vendor/

cp node_modules/topojson/build/topojson.js vendor/