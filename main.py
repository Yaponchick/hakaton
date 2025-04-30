# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from math import radians, cos, sin, asin, sqrt
from datetime import datetime, timezone
import requests
import geopandas as gpd
try:
    import mplcursors
    MPLCURSORS_AVAILABLE = True
except ImportError:
    MPLCURSORS_AVAILABLE = False
    print("*"*60)
    print("ПРЕДУПРЕЖДЕНИЕ: Библиотека 'mplcursors' не найдена.")
    print("Для интерактивных подсказок установите её: pip install mplcursors")
    print("Интерактивные подсказки будут отключены.")
    print("*"*60)

def haversine(lon1, lat1, lon2, lat2):
    """Рассчитывает расстояние между двумя точками на сфере"""
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # радиус Земли в км
    return c * r

def is_point_near_osm_stop(lon, lat, osm_stops, threshold_meters=50):
    for stop in osm_stops:
        distance_km = haversine(lon, lat, stop['lon'], stop['lat'])
        if distance_km * 1000 <= threshold_meters:
            return True
    return False

def load_route_data(excel_file):
    try:
        df = pd.read_excel(excel_file)
        route_data = {}
        for _, row in df.iterrows():
            if pd.notna(row['ID']):
                try:
                    route_id = int(row['ID'])
                    route_data[route_id] = {
                        'Название': row.get('Маршрут', 'Не указано'),
                        'Расстояние': row.get('Расстояние, км.', 'Не указано'),
                        'Время старта': row.get('Время старта замера, чч:мм:сек', 'Не указано'),
                        'Продолжительность': row.get('Продолжительность поездки, чч:мм:сек.', 'Не указано'),
                        'Средняя скорость': row.get('Средняя скорость, км/ч.', 'Не указано'),
                        'ФИО': row.get('Ф.И.О.', 'Не указано')
                    }
                except Exception:
                    pass
        return route_data
    except Exception as e:
        print(f"Ошибка при чтении Excel файла '{excel_file}': {e}")
        return {}

def parse_gpx(file_path):
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        namespaces = {'gpx': 'http://www.topografix.com/GPX/1/1'}
        points_data = []
        trkpts = root.findall('.//gpx:trk/gpx:trkseg/gpx:trkpt', namespaces)
        if not trkpts:
            trkpts = root.findall('.//gpx:trkpt', namespaces)
        for trkpt in trkpts:
            try:
                lat = float(trkpt.attrib['lat'])
                lon = float(trkpt.attrib['lon'])
                time_str = trkpt.find('gpx:time', namespaces)
                if time_str is None or not time_str.text:
                    continue
                dt_obj = datetime.fromisoformat(time_str.text.replace('Z', '+00:00'))
                if dt_obj.tzinfo is None:
                    dt_obj = dt_obj.replace(tzinfo=timezone.utc)
                points_data.append((lon, lat, dt_obj))
            except Exception:
                pass
        points_data.sort(key=lambda x: x[2])  # Сортировка по времени
        return points_data
    except Exception as e:
        print(f"Ошибка при парсинге GPX файла '{file_path}': {e}")
        return []

def get_osm_stops(min_lat, min_lon, max_lat, max_lon):
    overpass_url = "http://overpass-api.de/api/interpreter"
    query = f"""
    [out:json][timeout:60];
    (
      node["public_transport"="platform"]({min_lat},{min_lon},{max_lat},{max_lon});
      node["highway"="bus_stop"]({min_lat},{min_lon},{max_lat},{max_lon});
    );
    out center;
    """
    print(f"[OSM] Запрос остановок в области: S={min_lat:.4f}, W={min_lon:.4f}, N={max_lat:.4f}, E={max_lon:.4f}")
    stops = []
    try:
        response = requests.get(overpass_url, params={'data': query})
        response.raise_for_status()
        data = response.json()
        print(f"[OSM] Получено {len(data.get('elements', []))} элементов.")
        for element in data.get('elements', []):
            lat = element.get('lat')
            lon = element.get('lon')
            if lat is None or lon is None:
                center = element.get('center')
                if center:
                    lat = center.get('lat')
                    lon = center.get('lon')
            if lat is not None and lon is not None:
                stops.append({'lat': lat, 'lon': lon})
        print(f"[OSM] Найдено {len(stops)} валидных остановок.")
        return stops
    except requests.exceptions.RequestException as e:
        print(f"[OSM] Ошибка сети при запросе остановок: {e}")
    except Exception as e:
        print(f"[OSM] Неожиданная ошибка при обработке ответа OSM: {e}")
    return []

def plot_tracks_with_stops(all_tracks_data, route_data, uds_path=None):
    fig, ax = plt.subplots(figsize=(12, 10))
    plt.subplots_adjust(left=0.08, right=0.78, top=0.95, bottom=0.1)

    if uds_path:
        try:
            # Загрузка дорожной сети
            UDS = gpd.read_file(uds_path)
            print("SHP-файл успешно загружен!")
            # Перевод в EPSG:4326 (широта/долгота)
            UDS = UDS.to_crs(epsg=4326)
            # Отображение дорожной сети
            UDS.plot(ax=ax, edgecolor='grey', facecolor='none', linewidth=0.5, label='Дорожная сеть')
        except Exception as e:
            print(f"Ошибка при чтении SHP-файла: {e}")

    line_elements_map = {}
    visibility_state = {}
    has_data_to_plot = False
    all_points_exist = False
    min_lat_all, max_lat_all = 90.0, -90.0
    min_lon_all, max_lon_all = 180.0, -180.0
    global_show_osm_stops = True
    osm_stop_marker = None

    if 'tracks' in all_tracks_data:
        for track_info in all_tracks_data['tracks']:
            points_data = track_info['coords']
            if not points_data or len(points_data) < 2:
                continue
            longitudes = [p[0] for p in points_data]
            latitudes = [p[1] for p in points_data]
            min_lat_all = min(min_lat_all, min(latitudes))
            max_lat_all = max(max_lat_all, max(latitudes))
            min_lon_all = min(min_lon_all, min(longitudes))
            max_lon_all = max(max_lon_all, max(longitudes))
            all_points_exist = True

    osm_stops = []
    if all_points_exist:
        margin = 0.01
        osm_stops = get_osm_stops(
            min_lat_all - margin, min_lon_all - margin,
            max_lat_all + margin, max_lon_all + margin
        )

    unique_labels_seen = {}

    if 'tracks' in all_tracks_data:
        for idx, track_info in enumerate(all_tracks_data['tracks']):
            filename = track_info['filename']
            points_data = track_info['coords']
            base_filename = os.path.splitext(filename)[0]
            if not points_data or len(points_data) < 2:
                continue
            longitudes = [p[0] for p in points_data]
            latitudes = [p[1] for p in points_data]
            timestamps = [p[2] for p in points_data]

            total_distance_km = sum(
                haversine(longitudes[i - 1], latitudes[i - 1], longitudes[i], latitudes[i])
                for i in range(1, len(points_data))
            )
            duration_seconds = (timestamps[-1] - timestamps[0]).total_seconds()
            total_time_hours = duration_seconds / 3600.0 if duration_seconds > 0 else 0
            calculated_avg_speed_kmh = total_distance_km / total_time_hours if total_time_hours > 0 else 0.0

            try:
                route_id_str = base_filename.split('_')[0]
                route_id = int(route_id_str)
                route_info_excel = route_data.get(route_id, {})
                plot_label_base = f"{base_filename}"
            except ValueError:
                plot_label_base = base_filename
                route_info_excel = {}

            original_label = plot_label_base
            counter = 1
            while plot_label_base in unique_labels_seen:
                plot_label_base = f"{original_label}_{counter}"
                counter += 1
            unique_labels_seen[plot_label_base] = True

            line, = ax.plot(longitudes, latitudes, marker='.', linestyle='-', markersize=3,
                            label=plot_label_base, picker=5, linewidth=2)
            start_marker, = ax.plot(longitudes[0], latitudes[0], 'go', markersize=6, label='_nolegend_')
            end_marker, = ax.plot(longitudes[-1], latitudes[-1], 'ro', markersize=6, label='_nolegend_')

            valid_stop_points_lons = []
            valid_stop_points_lats = []
            confirmed_stop_points_lons = []
            confirmed_stop_points_lats = []
            unconfirmed_stop_points_lons = []
            unconfirmed_stop_points_lats = []

            speed_threshold_kmh = 1

            for i in range(1, len(points_data)):
                prev_lon, prev_lat, prev_time = points_data[i - 1]
                curr_lon, curr_lat, curr_time = points_data[i]
                delta_time_sec = (curr_time - prev_time).total_seconds()
                if delta_time_sec > 1e-6:
                    distance_km = haversine(prev_lon, prev_lat, curr_lon, curr_lat)
                    speed_kmh = (distance_km / delta_time_sec) * 3600.0
                    if speed_kmh < speed_threshold_kmh:
                        valid_stop_points_lons.append(curr_lon)
                        valid_stop_points_lats.append(curr_lat)
                        if is_point_near_osm_stop(curr_lon, curr_lat, osm_stops, threshold_meters=50):
                            confirmed_stop_points_lons.append(curr_lon)
                            confirmed_stop_points_lats.append(curr_lat)
                        else:
                            unconfirmed_stop_points_lons.append(curr_lon)
                            unconfirmed_stop_points_lats.append(curr_lat)

            stops_plot = None
            confirmed_plot = None
            unconfirmed_plot = None

            if valid_stop_points_lons:
                stops_plot, = ax.plot(valid_stop_points_lons, valid_stop_points_lats,
                                      'yo', markersize=8, linestyle='none', label='_nolegend_', picker=5)
            if confirmed_stop_points_lons:
                confirmed_plot, = ax.plot(confirmed_stop_points_lons, confirmed_stop_points_lats,
                                          'go', markersize=6, linestyle='none', label='_nolegend_', picker=5)
            if unconfirmed_stop_points_lats:
                unconfirmed_plot, = ax.plot(unconfirmed_stop_points_lons, unconfirmed_stop_points_lats,
                                            'ro', markersize=6, linestyle='none', label='_nolegend_', picker=5)

            line_elements_map[line] = {
                'label_base': plot_label_base,
                'start': start_marker,
                'end': end_marker,
                'stops': stops_plot,
                'confirmed_stops': confirmed_plot,
                'unconfirmed_stops': unconfirmed_plot,
                'calculated_speed': calculated_avg_speed_kmh,
                'calculated_distance': total_distance_km,
                'coords': points_data
            }
            visibility_state[plot_label_base] = True
            has_data_to_plot = True

    if not has_data_to_plot:
        print("Нет данных треков для отображения на графике.")
        plt.close(fig)
        return

    if osm_stops:
        osm_lon = [stop['lon'] for stop in osm_stops]
        osm_lat = [stop['lat'] for stop in osm_stops]
        osm_stop_marker, = ax.plot(osm_lon, osm_lat, 'rs', markersize=5, linestyle='none', label='ОСТАНОВКИ (OSM)', alpha=0.7)
        osm_stop_marker.set_visible(global_show_osm_stops)
    else:
        osm_stop_marker = None

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend = ax.legend(by_label.values(), by_label.keys(),
                       bbox_to_anchor=(1.03, 1), loc='upper left',
                       borderaxespad=0., prop={'size': 9})

    for leg_handle in legend.legend_handles:
        leg_handle.set_picker(True)
        leg_handle.set_pickradius(5)

    def on_pick(event):
        leg_artist = event.artist
        clicked_label = None
        legend_texts = legend.get_texts()
        labels = [text.get_text() for text in legend_texts]
        for i, text in enumerate(legend_texts):
            if leg_artist == legend.legend_handles[i]:
                clicked_label = text.get_text()
                break
        if clicked_label and clicked_label in visibility_state:
            new_visibility = not visibility_state[clicked_label]
            visibility_state[clicked_label] = new_visibility
            for line, data in line_elements_map.items():
                if data['label_base'] == clicked_label:
                    line.set_visible(new_visibility)
                    data['start'].set_visible(new_visibility)
                    data['end'].set_visible(new_visibility)
                    if data['stops']:
                        data['stops'].set_visible(new_visibility)
                    if data['confirmed_stops']:
                        data['confirmed_stops'].set_visible(new_visibility)
                    if data['unconfirmed_stops']:
                        data['unconfirmed_stops'].set_visible(new_visibility)
            leg_artist.set_alpha(1.0 if new_visibility else 0.2)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('pick_event', on_pick)

    ax_button_deselect = fig.add_axes([0.8, 0.015, 0.15, 0.04])
    button_deselect = Button(ax_button_deselect, 'Снять все')

    def deselect_all(event):
        for line, data in line_elements_map.items():
            base_label = data['label_base']
            line.set_visible(False)
            data['start'].set_visible(False)
            data['end'].set_visible(False)
            if data['stops']:
                data['stops'].set_visible(False)
            if data['confirmed_stops']:
                data['confirmed_stops'].set_visible(False)
            if data['unconfirmed_stops']:
                data['unconfirmed_stops'].set_visible(False)
            visibility_state[base_label] = False
        for handle in legend.legend_handles:
            handle.set_alpha(0.2)
        fig.canvas.draw_idle()

    button_deselect.on_clicked(deselect_all)

    ax_button_show_hide = fig.add_axes([0.64, 0.015, 0.15, 0.04])
    button_show_hide = Button(ax_button_show_hide, 'Скрыть остановки')

    def toggle_osm_stops(event):
        nonlocal global_show_osm_stops
        if osm_stop_marker is not None:
            global_show_osm_stops = not global_show_osm_stops
            osm_stop_marker.set_visible(global_show_osm_stops)
            button_show_hide.label.set_text('Скрыть остановки' if global_show_osm_stops else 'Показать остановки')
            fig.canvas.draw_idle()

    button_show_hide.on_clicked(toggle_osm_stops)

    if MPLCURSORS_AVAILABLE:
        cursor = mplcursors.cursor(list(line_elements_map.keys()), hover=True)

        @cursor.connect("add")
        def on_add(sel):
            origline = sel.artist
            if origline in line_elements_map:
                elements_data = line_elements_map[origline]
                label = elements_data['label_base']
                calc_speed = elements_data.get('calculated_speed', 0.0)
                calc_dist = elements_data.get('calculated_distance', 0.0)
                points_data = elements_data.get('coords', [])
                index = min(range(len(points_data)), key=lambda i: (points_data[i][0]-sel.target[0])**2 + (points_data[i][1]-sel.target[1])**2)
                if index > 0 and len(points_data) > index:
                    prev_point = points_data[index - 1]
                    curr_point = points_data[index]
                    delta_time_sec = (curr_point[2] - prev_point[2]).total_seconds()
                    if delta_time_sec > 0:
                        distance_km = haversine(prev_point[0], prev_point[1], curr_point[0], curr_point[1])
                        point_speed_kmh = (distance_km / delta_time_sec) * 3600
                    else:
                        point_speed_kmh = 0.0
                else:
                    point_speed_kmh = 0.0
                duration_seconds = (points_data[-1][2] - points_data[0][2]).total_seconds()
                hours, remainder = divmod(int(duration_seconds), 3600)
                minutes, seconds = divmod(remainder, 60)
                calc_dur_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
                tooltip_text = (
                    f"{label}\n"
                    f"Lon: {sel.target[0]:.5f}, Lat: {sel.target[1]:.5f}\n"
                    f"Расст.: {calc_dist:.2f} км\n"
                    f"Время: {calc_dur_str}\n"
                    f"Ср. скор.: {calc_speed:.1f} км/ч\n"
                    f"Скор. на точке: {point_speed_kmh:.1f} км/ч"
                )
                sel.annotation.set(text=tooltip_text)
                sel.annotation.get_bbox_patch().set(alpha=0.8)

    ax.set_title("Треки GPX с анализом остановок (кликните легенду для фильтра)")
    ax.set_xlabel("Долгота")
    ax.set_ylabel("Широта")
    ax.grid(True)
    plt.show()


if __name__ == "__main__":
    import os
    import geopandas as gpd  # Убедитесь, что библиотека импортирована

    gpx_directory = "tracks"
    excel_file = "Результаты треков ЛИМБ-21-1.xlsx"

    uds_path = 'Graph_Irkutsk_link/Graph_Irkutsk_link.SHP'

    if not os.path.isdir(gpx_directory):
        print(f"Ошибка: Директория '{gpx_directory}' не найдена.")
        exit()

    print(f"Загрузка данных из Excel: {excel_file}")
    route_data = load_route_data(excel_file)

    print(f"Поиск GPX файлов в директории: {gpx_directory}")
    gpx_files = [os.path.join(gpx_directory, f) for f in os.listdir(gpx_directory) if f.lower().endswith('.gpx')]

    if not gpx_files:
        print(f"В директории '{gpx_directory}' не найдено ни одного GPX файла (*.gpx).")
    else:
        print(f"Найдено {len(gpx_files)} GPX файлов. Начинаем обработку...")
        all_tracks_data = {'tracks': []}
        processed_files = 0
        for gpx_file in gpx_files:
            points_data = parse_gpx(gpx_file)
            if points_data:
                all_tracks_data['tracks'].append({
                    'filename': os.path.basename(gpx_file),
                    'coords': points_data
                })
                processed_files += 1
        print(
            f"Обработка файлов завершена. Успешно обработано (с точками времени): {processed_files} из {len(gpx_files)}.")

        if all_tracks_data['tracks']:
            plot_tracks_with_stops(all_tracks_data, route_data, uds_path=uds_path)
        else:
            print("Нет данных треков (с точками времени) для отображения после обработки файлов.")