<?php
$counterFile = 'counter.txt';

// Check if the counter file exists
if (!file_exists($counterFile)) {
    // Create the file and initialize it with 0
    file_put_contents($counterFile, 0);
}

// Read the current count
$count = (int)file_get_contents($counterFile);

// Increment the count
$count++;

// Write the new count back to the file
file_put_contents($counterFile, $count);

echo "This page has been visited $count times.";
?>
