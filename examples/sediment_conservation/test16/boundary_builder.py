#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 17:01:52 2026

@author: fgalifornisilv
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from typing import Iterable, List, Tuple

import pandas as pd


def build_synthetic_water_level_boundary(
	start_datetime: datetime | None = None,
	step_minutes: int = 10,
	pattern: Iterable[Tuple[str, int, float]] | None = None,
	swl_path: str | None = None,
	month_path: str | None = None,
) -> pd.DataFrame:
	"""
	Build a synthetic water-level boundary with equidistant 10-minute steps.

	Parameters
	----------
	start_datetime : datetime | None
		Start datetime for the boundary. If None, defaults to 2021-01-01 UTC.
	step_minutes : int
		Time step in minutes (default 10).
	pattern : Iterable[Tuple[str, int, float]] | None
		Sequence of (unit, count, level) entries. unit is "month" or "day".
		If None, uses the requested default pattern.
	swl_path : str | None
		If provided, saves datetime and swl to a .txt file (tab-separated).
	month_path : str | None
		If provided, saves datetime and month to a .txt file (tab-separated).

	Returns
	-------
	pd.DataFrame
		DataFrame with columns: datetime, swl, month.
	"""

	if start_datetime is None:
		start_datetime = datetime(2021, 1, 1, tzinfo=timezone.utc)
	elif start_datetime.tzinfo is None:
		start_datetime = start_datetime.replace(tzinfo=timezone.utc)

	if pattern is None:
		pattern = [
			("month", 1, 0.0),
			("day", 1, 5.0),
			("month", 1, 0.0),
			("day", 1, 7.0),
			("month", 1, 0.0),
			("day", 1, 10.0),
			("month", 15, 0.0),
		]

	if step_minutes <= 0:
		raise ValueError("step_minutes must be positive.")

	timestamps: List[pd.Timestamp] = []
	swl_values: List[float] = []

	current = pd.Timestamp(start_datetime)
	step = f"{int(step_minutes)}min"

	for unit, count, level in pattern:
		unit_lower = unit.lower()
		if unit_lower == "month":
			end = current + pd.DateOffset(months=int(count))
		elif unit_lower == "day":
			end = current + timedelta(days=int(count))
		else:
			raise ValueError(f"Unsupported unit '{unit}'. Use 'month' or 'day'.")

		segment_times = pd.date_range(
			start=current,
			end=end,
			freq=step,
			inclusive="left",
		)

		timestamps.extend(segment_times.to_list())
		swl_values.extend([float(level)] * len(segment_times))

		current = end

	df = pd.DataFrame(
		{
			"datetime": pd.to_datetime(timestamps, utc=True),
			"swl": swl_values,
		}
	)
	df["month"] = df["datetime"].dt.month

	if swl_path:
		df[["datetime", "swl"]].to_csv(swl_path, sep="\t", index=False)
	if month_path:
		df[["datetime", "month"]].to_csv(month_path, sep="\t", index=False)

	return df


if __name__ == "__main__":
	boundary = build_synthetic_water_level_boundary(
		start_datetime=datetime(2021, 1, 11, tzinfo=timezone.utc),
		swl_path="synthetic_swl.txt",
		month_path="synthetic_month.txt",
	)
	print(boundary.head())

	import matplotlib.pyplot as plt

	plt.figure(figsize=(12, 5))
	plt.plot(boundary["datetime"], boundary["swl"], linewidth=1.5)
	plt.xlabel("DateTime")
	plt.ylabel("Water Level (swl)")
	plt.title("Synthetic Water Level Boundary")
	plt.grid(True, alpha=0.3)
	plt.xticks(rotation=45)
	plt.tight_layout()
	plt.show()